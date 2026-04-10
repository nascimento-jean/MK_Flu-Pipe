#!/bin/bash
# =============================================================================
# MK Flu-Pipe — Montagem e Tipagem de Genomas de Influenza
# Versão: 1.6.2
# Desenvolvido por: Jean Phellipe Marques do Nascimento
# Laboratório de Vigilância Genômica — LACEN/AL
# =============================================================================
#
# NOVIDADES v1.6.2:
#   • Correção de subtipo_HA na Etapa 4 — agora extrai apenas o componente H
#     (ex: H5, H1, H3) em vez do subtipo completo HxNy
#   • get_nextclade_key() atualizada para aceitar subtipo_ha como Hx simples
#     (H1 → H1N1, H3 → H3N2, H5 → H5 all-clades)
#   • Suporte a H5Nx no Nextclade — dataset comunitário:
#       community/moncla-lab/iav-h5/ha/all-clades
#     Cobre todos os clados H5 (H5N1, H5N8, H5N6…); exibe apenas "clade"
#     ex: 2.3.4.4b | 2.2.1
#   • Nova coluna classificacao_BLAST no typing_results.tsv (após subtipo_NA):
#     concatena subtipo_HA + subtipo_NA → ex: H5N8, H1N1, H3N2
#     Para Flu B exibe a linhagem (Victoria/Yamagata)
#   • Nova coluna segmentos no typing_results.tsv (após hit_blast_NA):
#     lista os segmentos presentes no multifasta de assembly_final/
#     separados por "|" → ex: 1|2|3|4|5|6|7|8  ou  2|4|8
#       nextstrain/flu/b/ha/KX058884  (classifica Victoria e Yamagata juntos)
#     Os datasets separados B_Victoria e B_Yamagata foram removidos.
#   • clado_nextclade para FluB exibe "clade/lineage"
#     ex: V1A.3a.2/Victoria  |  Y3/Yamagata
#   • get_nextclade_key() simplificada: qualquer FluB → chave "B"
#   • parse_nextclade_clade() detecta automaticamente se o TSV contém
#     "lineage" (FluB) ou "short-clade" (FluA) e usa a coluna correta
#   • Suporte multi-plataforma de sequenciamento:
#       - Illumina paired-end padrão BaseSpace (*_L001_R1_001.fastq.gz)
#       - Paired-end genérico (*_R1*.fastq.gz, *_1.fastq.gz, *_R1_*.fastq.gz)
#       - Single-end (ONT, Ion Torrent, PacBio, e outros)
#       - SRA/NCBI: paired (*_1.fastq.gz + *_2.fastq.gz) e single (*.fastq.gz)
#     Detecção automática do modo de entrada; opção de forçar via --mode
#
# NOVIDADES v1.6.1:
#   • Correção na Etapa 5 — cabeçalho do TSV de fallback agora inclui a
#     coluna "short-clade", compatível com o TSV real gerado pelo Nextclade
#   • Correção na Etapa 6 — clado_nextclade agora lê "clade" e "short-clade"
#     pelo nome da coluna (não por posição) e exibe no formato "clade/short-clade"
#     ex: J.2.3/2a.3a.1  |  5a.2a/6B.1A.5a.2a
#   • Para subtipos sem dataset Nextclade (H5N1, H7N9 etc.) as colunas
#     clado_nextclade e qc_nextclade ficam "—" e fonte_classificacao = "BLAST"
#
# NOVIDADES v1.6.0:
#   • Substituição do ABRicate pelo pipeline BLAST + Nextclade:
#       Etapa 3 — Banco NCBI Influenza (download automático / atualização
#                 trimestral em ~/mk_flupipe_db/)
#       Etapa 4 — BLASTN nos segmentos HA (seg.4) e NA (seg.6) para
#                 determinar tipo (A/B/C/D) e subtipo (HxNy)
#       Etapa 5 — Nextclade CLI para classificação de clado nos subtipos
#                 H1N1pdm, H3N2 e Influenza B (Victoria + Yamagata)
#       Etapa 6 — Consolidação em typing_results.tsv por corrida
#
# NOVIDADES v1.5.0:
#   • Sistema de checkpoint/resume — retoma de onde parou após interrupção.
#     Arquivos de marcação em $output_dir/.checkpoints/
#
# Ambiente conda : mk_flu (blast + nextclade)
# Banco          : ~/mk_flupipe_db/influenza.fna (NCBI Influenza DB)
# Atualização    : automática a cada 90 dias
# =============================================================================

set -euo pipefail

# --------------------------------------------------------------------------- #
# Configurações globais
# --------------------------------------------------------------------------- #

DB_DIR="$HOME/mk_flupipe_db"
DB_FILE="$DB_DIR/influenza.fna"
DB_BLAST="$DB_DIR/influenza_blast_db"
DB_TIMESTAMP="$DB_DIR/.last_update"
DB_MAX_DAYS=90
DB_URL="https://ftp.ncbi.nih.gov/genomes/INFLUENZA/influenza.fna.gz"

CONDA_ENV_MKFLU="mk_flu"

# Modo de sequenciamento (detectado automaticamente ou forçado via --mode)
# Valores: "illumina_paired" | "generic_paired" | "single" | "sra_paired" | "sra_single"
SEQ_MODE=""

# Datasets Nextclade disponíveis para influenza humana
# Chave = identificador interno; valor = nome do dataset no Nextclade
# FluB: dataset unificado que classifica Victoria e Yamagata juntos,
#        usando as colunas "clade" e "lineage" no TSV de saída.
declare -A NEXTCLADE_DATASETS=(
    ["H1N1"]="nextstrain/flu/h1n1pdm/ha/MW626062"
    ["H3N2"]="nextstrain/flu/h3n2/ha/EPI1857216"
    ["B"]="nextstrain/flu/b/ha/KX058884"
    ["H5"]="community/moncla-lab/iav-h5/ha/all-clades"
)

# --------------------------------------------------------------------------- #
# Funções auxiliares
# --------------------------------------------------------------------------- #

timestamp() { date '+%Y-%m-%d %H:%M:%S'; }

log_info()    { echo "[$(timestamp)] [INFO]    $*"; }
log_warn()    { echo "[$(timestamp)] [AVISO]   $*"; }
log_error()   { echo "[$(timestamp)] [ERRO]    $*" >&2; }
log_skip()    { echo "[$(timestamp)] [SKIP]    $*"; }
log_section() {
    echo ""
    echo "════════════════════════════════════════════════════════════════════════"
    echo "  $*"
    echo "════════════════════════════════════════════════════════════════════════"
}

check_dependency() {
    local cmd="$1"
    if ! command -v "$cmd" &>/dev/null; then
        log_error "Dependência não encontrada: '$cmd'."
        exit 1
    fi
}

mark_done() { touch "$1"; }
is_done()   { [ -f "$1" ]; }

find_conda_base() {
    for candidate in \
        "$HOME/miniconda3" "$HOME/miniforge3" "$HOME/mambaforge" \
        "$HOME/anaconda3"  "/opt/conda"       "/opt/miniconda3"  \
        "/opt/miniforge3"
    do
        [ -f "$candidate/bin/activate" ] && echo "$candidate" && return 0
    done
    if [ -n "${CONDA_EXE:-}" ]; then
        dirname "$(dirname "$CONDA_EXE")" && return 0
    fi
    return 1
}

days_since_modified() {
    local file="$1"
    local now mod
    now=$(date +%s)
    mod=$(date -r "$file" +%s 2>/dev/null || echo 0)
    echo $(( (now - mod) / 86400 ))
}

# Extrai tipo (A/B/C/D) e subtipo do título do melhor hit BLAST do NCBI Influenza DB
# O header NCBI Influenza tem formato variado, ex:
#   Influenza A virus (A/California/07/2009(H1N1)) segment 4 ...
#   gb|CY121680|Influenza B virus (B/Brisbane/60/2008) segment 4 ...
parse_blast_hit() {
    local stitle="$1"
    local tipo subtipo

    if echo "$stitle" | grep -qiE "influenza A|type A|\(A/"; then
        tipo="A"
    elif echo "$stitle" | grep -qiE "influenza B|type B|\(B/"; then
        tipo="B"
    elif echo "$stitle" | grep -qiE "influenza C|\(C/"; then
        tipo="C"
    elif echo "$stitle" | grep -qiE "influenza D|\(D/"; then
        tipo="D"
    else
        tipo="desconhecido"
    fi

    # Tenta capturar HxNy
    subtipo=$(echo "$stitle" | grep -oiE 'H[0-9]{1,2}N[0-9]{1,2}' | head -1 || true)

    # Para Flu B sem HxNy, verifica linhagem
    if [ "$tipo" = "B" ] && [ -z "$subtipo" ]; then
        if echo "$stitle" | grep -qi "Victoria\|Vic"; then
            subtipo="Victoria"
        elif echo "$stitle" | grep -qi "Yamagata\|Yam"; then
            subtipo="Yamagata"
        else
            subtipo="nd"
        fi
    fi

    echo "${tipo}|${subtipo:-nd}"
}

# Mapeia tipo + subtipo_HA para a chave do dataset Nextclade
# FluA: H1 → "H1N1" | H3 → "H3N2" | H5 (qualquer Nx) → "H5" | outros → ""
# FluB: qualquer linhagem → "B" (dataset unificado Victoria + Yamagata)
# H5:   dataset comunitário cobre todos os clados H5Nx (H5N1, H5N8, H5N6…)
get_nextclade_key() {
    local tipo="$1"
    local subtipo_ha="$2"
    local sub_upper="${subtipo_ha^^}"

    if [ "$tipo" = "A" ]; then
        case "$sub_upper" in
            H1N1*|H1) echo "H1N1" ;;
            H3N2*|H3) echo "H3N2" ;;
            H5*)       echo "H5"   ;;
            *)         echo "" ;;
        esac
    elif [ "$tipo" = "B" ]; then
        # Dataset unificado cobre Victoria e Yamagata — chave sempre "B"
        echo "B"
    else
        echo ""
    fi
}

# =============================================================================
# FUNÇÃO: extrai clado a partir do TSV Nextclade, adaptando-se ao tipo viral
#
# FluA (H1N1, H3N2) → lê "clade" e "short-clade" → formato "clade/short-clade"
#                       ex: J.2.3/2a.3a.1 | 3C.2a1b/2a
# FluB (dataset unificado) → lê "clade" e "lineage" → formato "clade/lineage"
#                       ex: V1A.3a.2/Victoria | Y3/Yamagata
# H5Nx (all-clades)    → lê apenas "clade" → formato "clade"
#                       ex: 2.3.4.4b | 2.2.1
#
# A detecção é automática pelo cabeçalho do TSV:
#   presença de "lineage"     → FluB
#   presença de "short-clade" → FluA (H1N1/H3N2)
#   apenas "clade"            → H5Nx (ou qualquer dataset sem second field)
#
# Retorna "—" quando não há classificação válida (subtipo atípico,
# arquivo ausente, ou falha do Nextclade).
# =============================================================================
parse_nextclade_clade() {
    local nc_file="$1"

    awk -F'\t' '
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                if ($i == "clade")       col_clade = i
                if ($i == "short-clade") col_short = i
                if ($i == "lineage")     col_lineage = i
            }
            next
        }
        NR == 2 {
            clade = (col_clade   ? $col_clade   : "")
            short = (col_short   ? $col_short   : "")
            lin   = (col_lineage ? $col_lineage : "")

            # Valores inválidos gerados pelo fallback interno
            if (clade == "N/A" || clade == "ERRO" || clade == "") {
                print "\xe2\x80\x94"   # "—"
                exit
            }

            # FluB: usa "lineage" como segundo campo
            if (col_lineage && lin != "" && lin != "N/A") {
                print clade "/" lin
                exit
            }

            # FluA: usa "short-clade" como segundo campo
            if (col_short && short != "" && short != "N/A" && short != clade) {
                print clade "/" short
                exit
            }

            # Fallback: apenas o clade sem segundo campo
            print clade
        }
    ' "$nc_file" 2>/dev/null || echo "—"
}

# --------------------------------------------------------------------------- #
# Detecção automática do modo de sequenciamento / entrada
# --------------------------------------------------------------------------- #
#
# Prioridade de detecção (da mais específica para a mais genérica):
#   1. illumina_paired  — *_L001_R1_001.fastq.gz  (BaseSpace padrão)
#   2. sra_paired       — *_1.fastq.gz + *_2.fastq.gz  (SRA paired)
#   3. generic_paired   — *_R1*.fastq.gz | *_R1_*.fastq.gz  (paired genérico)
#   4. single           — qualquer *.fastq.gz restante (single-end)
#   5. sra_single       — *.fastq.gz único sem par identificável
#
# O modo pode ser forçado pelo usuário com --mode <modo>.
# --------------------------------------------------------------------------- #

detect_seq_mode() {
    local dir="$1"
    local forced="${2:-}"

    if [ -n "$forced" ]; then
        SEQ_MODE="$forced"
        log_info "Modo de sequenciamento forçado pelo usuário: $SEQ_MODE"
        return 0
    fi

    shopt -s nullglob

    # 1. Illumina BaseSpace padrão
    local ill=("$dir"/*_L001_R1_001.fastq.gz)
    if [ ${#ill[@]} -gt 0 ]; then
        SEQ_MODE="illumina_paired"
        shopt -u nullglob
        return 0
    fi

    # 2. SRA paired: arquivos terminando em _1.fastq.gz com correspondente _2
    local sra1=("$dir"/*_1.fastq.gz)
    if [ ${#sra1[@]} -gt 0 ]; then
        local has_pair=0
        for f in "${sra1[@]}"; do
            base="${f%_1.fastq.gz}"
            [ -f "${base}_2.fastq.gz" ] && has_pair=1 && break
        done
        if [ "$has_pair" -eq 1 ]; then
            SEQ_MODE="sra_paired"
            shopt -u nullglob
            return 0
        fi
    fi

    # 3. Paired genérico: _R1_ ou _R1. (sem sufixo Illumina)
    local gen_r1=("$dir"/*_R1_*.fastq.gz "$dir"/*_R1.fastq.gz)
    if [ ${#gen_r1[@]} -gt 0 ]; then
        SEQ_MODE="generic_paired"
        shopt -u nullglob
        return 0
    fi

    # 4. Single-end genérico / SRA single (qualquer .fastq.gz restante)
    local any=("$dir"/*.fastq.gz)
    if [ ${#any[@]} -gt 0 ]; then
        SEQ_MODE="single"
        shopt -u nullglob
        return 0
    fi

    shopt -u nullglob
    SEQ_MODE=""
}

# Retorna (via stdout) os arquivos R1/único de cada amostra,
# um por linha, no formato: NOME_AMOSTRA<TAB>R1_FILE[<TAB>R2_FILE]
list_samples() {
    local dir="$1"
    shopt -s nullglob

    case "$SEQ_MODE" in

        illumina_paired)
            for r1 in "$dir"/*_L001_R1_001.fastq.gz; do
                local name
                name=$(basename "$r1" _L001_R1_001.fastq.gz)
                local r2="$dir/${name}_L001_R2_001.fastq.gz"
                echo -e "${name}\t${r1}\t${r2}"
            done
            ;;

        sra_paired)
            for r1 in "$dir"/*_1.fastq.gz; do
                local name
                name=$(basename "$r1" _1.fastq.gz)
                local r2="$dir/${name}_2.fastq.gz"
                echo -e "${name}\t${r1}\t${r2}"
            done
            ;;

        generic_paired)
            # Suporta _R1_ e _R1. como separadores
            for r1 in "$dir"/*_R1_*.fastq.gz "$dir"/*_R1.fastq.gz; do
                [ -f "$r1" ] || continue
                local name r2
                # Descobre o R2 correspondente
                r2="${r1/_R1_/_R2_}"
                r2="${r2/_R1.fastq.gz/_R2.fastq.gz}"
                name=$(basename "$r1")
                # Remove sufixo variável para obter nome da amostra
                name="${name/_R1*/}"
                name="${name%.fastq.gz}"
                echo -e "${name}\t${r1}\t${r2}"
            done
            ;;

        single|sra_single)
            for f in "$dir"/*.fastq.gz; do
                local name
                name=$(basename "$f" .fastq.gz)
                # Remove sufixos comuns de SRA single (_pass, _fail, etc.)
                name="${name%_pass}"
                name="${name%_fail}"
                echo -e "${name}\t${f}"
            done
            ;;
    esac

    shopt -u nullglob
}

# --------------------------------------------------------------------------- #
# Cabeçalho
# --------------------------------------------------------------------------- #

echo ""
echo "████████████████████████████████████████████████████████████████████████"
echo "█          SCRIPT PARA MONTAGEM DE GENOMAS DE INFLUENZA                █"
echo "█                        MK Flu-Pipe v1.6.2                            █"
echo "█              Desenvolvido por Jean Phellipe M. do Nascimento         █"
echo "█        LABORATORIO CENTRAL DE SAUDE PUBLICA DE ALAGOAS - LACEN/AL    █"
echo "████████████████████████████████████████████████████████████████████████"
echo ""
log_info "Início: $(timestamp)"
echo ""

# --------------------------------------------------------------------------- #
# Leitura de parâmetros
# --------------------------------------------------------------------------- #

input_dir="${1:-}"
output_dir="${2:-}"
irma_module="${3:-}"
forced_mode="${4:-}"   # opcional: illumina_paired | sra_paired | generic_paired | single

# Suporte a --mode flag
if [[ "$input_dir" == "--mode" ]]; then
    forced_mode="$output_dir"
    input_dir="${irma_module:-}"
    output_dir="${4:-}"
    irma_module="${5:-}"
fi

if [ -z "$input_dir" ] || [ -z "$output_dir" ] || [ -z "$irma_module" ]; then
    echo ">>> Modo interativo <<<"
    read -e -p "Caminho do diretório de entrada (FASTQ): " input_dir
    input_dir=$(echo "$input_dir" | xargs)
    read -e -p "Caminho do diretório de saída: " output_dir
    output_dir=$(echo "$output_dir" | xargs)
    echo "Módulos disponíveis: FLU | FLU-utr | FLU-lowQC | FLU_AD | FLU-minion"
    read -p "Módulo IRMA: " irma_module
    echo ""
    echo "Modos disponíveis (deixe em branco para detecção automática):"
    echo "  illumina_paired — Illumina BaseSpace (*_L001_R1_001.fastq.gz)"
    echo "  sra_paired      — SRA paired-end (*_1.fastq.gz + *_2.fastq.gz)"
    echo "  generic_paired  — Paired-end genérico (*_R1_*.fastq.gz)"
    echo "  single          — Single-end / ONT / Ion Torrent / PacBio"
    read -p "Modo de sequenciamento [auto]: " forced_mode
fi

# --------------------------------------------------------------------------- #
# Validações iniciais
# --------------------------------------------------------------------------- #

log_section "Validando parâmetros e dependências"

[ ! -d "$input_dir" ] && { log_error "Diretório de entrada não encontrado: '$input_dir'"; exit 1; }
mkdir -p "$output_dir" || { log_error "Não foi possível criar '$output_dir'"; exit 1; }

# Detecta (ou valida) o modo de sequenciamento
detect_seq_mode "$input_dir" "$forced_mode"

if [ -z "$SEQ_MODE" ]; then
    log_error "Nenhum arquivo .fastq.gz encontrado em '$input_dir'."
    log_error "Verifique se os arquivos estão no diretório correto."
    exit 1
fi

log_info "Entrada : $input_dir"
log_info "Saída   : $output_dir"
log_info "Módulo  : $irma_module"
log_info "Modo    : $SEQ_MODE"

check_dependency "IRMA"
check_dependency "free"
check_dependency "nproc"
check_dependency "wget"

# Valida ambiente conda mk_flu
conda_base=$(find_conda_base) || {
    log_error "Instalação do Conda não encontrada."
    exit 1
}
log_info "Conda base detectada: $conda_base"

if ! "$conda_base/bin/conda" env list | grep -qE "^${CONDA_ENV_MKFLU}[[:space:]]"; then
    log_error "Ambiente conda '${CONDA_ENV_MKFLU}' não encontrado."
    log_error "Crie com: conda create -n mk_flu -c conda-forge -c bioconda blast nextclade -y"
    exit 1
fi

CONDA_RUN="$conda_base/bin/conda run --no-capture-output -n $CONDA_ENV_MKFLU"

$CONDA_RUN blastn -version &>/dev/null || {
    log_error "BLAST não encontrado no ambiente '$CONDA_ENV_MKFLU'."; exit 1; }
$CONDA_RUN nextclade --version &>/dev/null || {
    log_error "Nextclade não encontrado no ambiente '$CONDA_ENV_MKFLU'."; exit 1; }

log_info "BLAST e Nextclade verificados no ambiente '${CONDA_ENV_MKFLU}'."

# --------------------------------------------------------------------------- #
# Checkpoint setup
# --------------------------------------------------------------------------- #

seg_base="$output_dir/assembly_final/segments"
checkpoint_dir="$output_dir/.checkpoints"
mkdir -p "$checkpoint_dir"

amostras_retomadas=0
shopt -s nullglob
done_files=("$checkpoint_dir"/*.irma.done)
shopt -u nullglob
amostras_retomadas=${#done_files[@]}

if [ "$amostras_retomadas" -gt 0 ]; then
    echo ""
    log_info "╔══════════════════════════════════════════════════════════════╗"
    log_info "║  MODO RETOMADA: $amostras_retomadas amostra(s) já processada(s).          ║"
    log_info "║  O pipeline continuará de onde parou.                        ║"
    log_info "╚══════════════════════════════════════════════════════════════╝"
    echo ""
fi

# --------------------------------------------------------------------------- #
# Contagem de amostras
# --------------------------------------------------------------------------- #

# Constrói lista de amostras em arquivo temporário para reutilização
samples_list_file="$output_dir/.samples_list.tsv"
mkdir -p "$output_dir"
list_samples "$input_dir" > "$samples_list_file"

total_amostras=$(wc -l < "$samples_list_file" | tr -d ' ')

if [ "$total_amostras" -eq 0 ]; then
    log_warn "Nenhuma amostra detectada no modo '$SEQ_MODE' em '$input_dir'."
    log_warn "Verifique se os arquivos .fastq.gz seguem o padrão esperado para o modo selecionado."
    exit 1
fi
log_info "Amostras encontradas: $total_amostras (modo: $SEQ_MODE)"

# --------------------------------------------------------------------------- #
# Etapa 1 — Montagem com IRMA
# --------------------------------------------------------------------------- #

log_section "Etapa 1 — Montagem de genomas com IRMA"

mkdir -p "$output_dir/assembly_final"
run_log_file="$output_dir/run_log.txt"

if [ "$amostras_retomadas" -gt 0 ]; then
    { echo ""; echo "=== RETOMADA em $(timestamp) ==="; } >> "$run_log_file"
else
    > "$run_log_file"
fi

amostra_atual=0; amostras_ok=0; amostras_puladas=0; amostras_falha=0

while IFS=$'\t' read -r file_name input_file_r1 input_file_r2; do
    checkpoint_irma="$checkpoint_dir/${file_name}.irma.done"
    amostra_atual=$((amostra_atual + 1))

    if is_done "$checkpoint_irma"; then
        log_skip "[$amostra_atual/$total_amostras] $file_name — IRMA já concluído. Pulando."
        amostras_ok=$((amostras_ok + 1))
        amostras_puladas=$((amostras_puladas + 1))
        continue
    fi

    log_info "[$amostra_atual/$total_amostras] Processando: $file_name"

    sample_output_dir="$output_dir/$file_name"
    mem_before=$(free -m | awk '/Mem/ {print $3}')
    mem_pct=$(free -m | awk '/Mem/ {printf "%.1f", $3/$2*100}')
    cores=$(nproc --all)

    {
        echo ""
        echo "=== Amostra: $file_name ==="
        echo "  Modo: $SEQ_MODE | Memória antes : ${mem_before} MB (${mem_pct}%) | Núcleos: $cores"
    } >> "$run_log_file"

    # ── Execução do IRMA conforme o modo ──────────────────────────────────
    irma_ok=0

    case "$SEQ_MODE" in

        illumina_paired|sra_paired|generic_paired)
            # Paired-end: verifica R2 antes de chamar IRMA
            if [ -z "$input_file_r2" ] || [ ! -f "$input_file_r2" ]; then
                log_warn "  R2 não encontrado para '$file_name'. Pulando."
                echo "[$file_name] AVISO: R2 ausente." >> "$run_log_file"
                amostras_falha=$((amostras_falha + 1))
                continue
            fi
            if IRMA "$irma_module" "$input_file_r1" "$input_file_r2" "$sample_output_dir" \
                    2>&1 | tee -a "$run_log_file"; then
                irma_ok=1
            fi
            ;;

        single|sra_single)
            # Single-end: apenas um arquivo de entrada
            if IRMA "$irma_module" "$input_file_r1" "$sample_output_dir" \
                    2>&1 | tee -a "$run_log_file"; then
                irma_ok=1
            fi
            ;;
    esac

    if [ "$irma_ok" -eq 1 ]; then
        amended_dir="$sample_output_dir/amended_consensus"
        if [ -d "$amended_dir" ] && compgen -G "$amended_dir/*.fa" > /dev/null; then
            cat "$amended_dir"/*.fa > "$output_dir/assembly_final/$file_name.fasta"
            log_info "  → Consenso: assembly_final/$file_name.fasta"
            mark_done "$checkpoint_irma"
            amostras_ok=$((amostras_ok + 1))
        else
            log_warn "  → amended_consensus vazio/ausente para '$file_name'."
            echo "  AVISO: amended_consensus ausente." >> "$run_log_file"
            amostras_falha=$((amostras_falha + 1))
        fi
    else
        log_error "  → IRMA falhou para '$file_name'."
        echo "  ERRO: IRMA retornou código != 0." >> "$run_log_file"
        amostras_falha=$((amostras_falha + 1))
    fi

    mem_after=$(free -m | awk '/Mem/ {print $3}')
    mem_pct2=$(free -m | awk '/Mem/ {printf "%.1f", $3/$2*100}')
    echo "  Memória depois: ${mem_after} MB (${mem_pct2}%)" >> "$run_log_file"
done < "$samples_list_file"

echo ""
log_info "Montagem — OK: $amostras_ok (${amostras_puladas} retomadas) | Falha: $amostras_falha"

if [ "$amostras_ok" -eq 0 ]; then
    log_error "Nenhuma amostra montada com sucesso. Abortando."
    exit 1
fi

# --------------------------------------------------------------------------- #
# Etapa 2 — Extração dos segmentos (1–8)
# --------------------------------------------------------------------------- #

log_section "Etapa 2 — Extração dos segmentos (1–8)"

checkpoint_segments="$checkpoint_dir/.segments.done"

if is_done "$checkpoint_segments"; then
    log_skip "Segmentos já extraídos (checkpoint). Pulando etapa 2."
else
    for seg in 1 2 3 4 5 6 7 8; do
        mkdir -p "$seg_base/single_segment${seg}"
        > "$seg_base/segmento_${seg}.fasta"
    done

    shopt -s nullglob
    fasta_files=("$output_dir/assembly_final"/*.fasta)
    shopt -u nullglob

    if [ ${#fasta_files[@]} -eq 0 ]; then
        log_warn "Nenhum .fasta em assembly_final. Pulando extração."
    else
        for arquivo in "${fasta_files[@]}"; do
            nome_amostra=$(basename "$arquivo" .fasta)
            for seg in 1 2 3 4 5 6 7 8; do
                grep -A 1 "_${seg}$" "$arquivo" | grep -v '^--$' \
                    >> "$seg_base/segmento_${seg}.fasta" || true
                grep -A 1 "_${seg}$" "$arquivo" | grep -v '^--$' \
                    > "$seg_base/single_segment${seg}/${nome_amostra}_segmento_${seg}.fasta" || true
            done
        done
        log_info "Segmentos extraídos para: $seg_base"
    fi

    # Substitui bases inválidas por N em todos os 8 segmentos
    for seg in 1 2 3 4 5 6 7 8; do
        seg_file="$seg_base/segmento_${seg}.fasta"
        [ -s "$seg_file" ] && sed -i '/^>/!s/[^ACTGactg]/N/g' "$seg_file"
        find "$seg_base/single_segment${seg}/" -type f -name "*.fasta" -size +0c \
            -exec sed -i '/^>/!s/[^ACTGactg]/N/g' {} \;
    done

    log_info "Bases inválidas → 'N' concluído (todos os 8 segmentos)."
    mark_done "$checkpoint_segments"
fi

# --------------------------------------------------------------------------- #
# Etapa 3 — Banco NCBI Influenza (download / atualização trimestral)
# --------------------------------------------------------------------------- #

log_section "Etapa 3 — Banco NCBI Influenza (~/mk_flupipe_db/)"

checkpoint_db="$checkpoint_dir/.ncbi_db.done"
mkdir -p "$DB_DIR"

_download_ncbi_db() {
    log_info "Baixando banco NCBI Influenza..."
    local gz_file="$DB_DIR/influenza.fna.gz"

    if ! wget -q --show-progress -O "$gz_file" "$DB_URL"; then
        log_warn "Falha no download. Sem conexão?"
        rm -f "$gz_file"
        return 1
    fi

    log_info "Descompactando..."
    gunzip -f "$gz_file"

    log_info "Formatando banco BLAST (makeblastdb)..."
    $CONDA_RUN makeblastdb \
        -in "$DB_FILE" \
        -dbtype nucl \
        -out "$DB_BLAST" \
        -title "NCBI_Influenza" \
        -parse_seqids \
        >> "$run_log_file" 2>&1

    date '+%Y-%m-%d' > "$DB_TIMESTAMP"
    log_info "Banco pronto em: $DB_DIR"
    return 0
}

if is_done "$checkpoint_db" && [ -f "$DB_FILE" ] && [ -f "${DB_BLAST}.nhr" ]; then
    log_skip "Banco NCBI já preparado nesta corrida (checkpoint). Pulando etapa 3."
else
    if [ ! -f "$DB_FILE" ] || [ ! -f "${DB_BLAST}.nhr" ]; then
        log_info "Banco não encontrado. Realizando download inicial..."
        _download_ncbi_db || {
            log_error "Não foi possível baixar o banco NCBI Influenza. Verifique a conexão."
            exit 1
        }
    else
        days_old=0
        [ -f "$DB_TIMESTAMP" ] && days_old=$(days_since_modified "$DB_TIMESTAMP") || days_old=$((DB_MAX_DAYS + 1))

        if [ "$days_old" -gt "$DB_MAX_DAYS" ]; then
            log_warn "Banco com ${days_old} dias (limite: ${DB_MAX_DAYS}). Tentando atualização..."
            if _download_ncbi_db; then
                log_info "Banco atualizado com sucesso."
            else
                log_warn "Atualização falhou. Continuando com banco existente (${days_old} dias)."
            fi
        else
            log_info "Banco NCBI OK — atualizado há ${days_old} dia(s)."
        fi
    fi
    mark_done "$checkpoint_db"
fi

# --------------------------------------------------------------------------- #
# Etapa 4 — BLASTN nos segmentos HA (seg.4) e NA (seg.6)
# --------------------------------------------------------------------------- #

log_section "Etapa 4 — BLASTN: tipagem por HA (seg.4) e NA (seg.6)"

checkpoint_blast="$checkpoint_dir/.blast.done"
blast_dir="$output_dir/assembly_final/blast_results"
mkdir -p "$blast_dir"
blast_summary="$blast_dir/blast_typing_summary.tsv"

BLAST_OUTFMT="6 qseqid sseqid pident length evalue bitscore stitle"
BLAST_PARAMS="-max_target_seqs 1 -max_hsps 1 -perc_identity 80 -evalue 1e-10"

if is_done "$checkpoint_blast"; then
    log_skip "BLAST já concluído (checkpoint). Pulando etapa 4."
else
    echo -e "amostra\ttipo_blast\tsubtipo_HA\tsubtipo_NA\thit_HA\thit_NA" > "$blast_summary"

    shopt -s nullglob
    ha_files=("$seg_base/single_segment4"/*.fasta)
    shopt -u nullglob

    for ha_file in "${ha_files[@]}"; do
        sample_name=$(basename "$ha_file" _segmento_4.fasta)

        # Retomada parcial: pula se já está no resumo
        if grep -q "^${sample_name}"$'\t' "$blast_summary" 2>/dev/null; then
            log_skip "  BLAST já realizado para $sample_name."
            continue
        fi

        log_info "  Processando: $sample_name"

        # ── BLAST HA ──────────────────────────────────────────────────────
        ha_blast_out="$blast_dir/${sample_name}_HA_blast.tsv"
        tipo_final="desconhecido"; subtipo_ha="nd"; hit_ha="sem_resultado"

        if [ -s "$ha_file" ]; then
            $CONDA_RUN blastn \
                -query "$ha_file" \
                -db "$DB_BLAST" \
                -out "$ha_blast_out" \
                -outfmt "$BLAST_OUTFMT" \
                $BLAST_PARAMS \
                -num_threads "$(nproc --all)" \
                >> "$run_log_file" 2>&1 || true

            if [ -s "$ha_blast_out" ]; then
                # stitle = todas as colunas a partir da 7ª
                best_stitle=$(awk 'NR==1{for(i=7;i<=NF;i++) printf "%s%s",$i,(i<NF?" ":"\n")}' "$ha_blast_out")
                parsed=$(parse_blast_hit "$best_stitle")
                tipo_final="${parsed%%|*}"
                subtipo_ha_raw="${parsed##*|}"
                # Extrai apenas Hx se vier embutido em HxNy
                subtipo_ha=$(echo "$subtipo_ha_raw" | grep -oiE 'H[0-9]{1,2}' | head -1 || echo "$subtipo_ha_raw")
                hit_ha=$(awk 'NR==1{print $2}' "$ha_blast_out")
                log_info "    HA → Tipo: $tipo_final | Subtipo HA: $subtipo_ha | Hit: $hit_ha"
            else
                log_warn "    HA → Sem hit BLAST para $sample_name."
            fi
        else
            log_warn "    HA → Arquivo vazio/ausente para $sample_name."
        fi

        # ── BLAST NA ──────────────────────────────────────────────────────
        na_file="$seg_base/single_segment6/${sample_name}_segmento_6.fasta"
        na_blast_out="$blast_dir/${sample_name}_NA_blast.tsv"
        subtipo_na="nd"; hit_na="sem_resultado"

        if [ -s "$na_file" ]; then
            $CONDA_RUN blastn \
                -query "$na_file" \
                -db "$DB_BLAST" \
                -out "$na_blast_out" \
                -outfmt "$BLAST_OUTFMT" \
                $BLAST_PARAMS \
                -num_threads "$(nproc --all)" \
                >> "$run_log_file" 2>&1 || true

            if [ -s "$na_blast_out" ]; then
                best_stitle_na=$(awk 'NR==1{for(i=7;i<=NF;i++) printf "%s%s",$i,(i<NF?" ":"\n")}' "$na_blast_out")
                parsed_na=$(parse_blast_hit "$best_stitle_na")
                subtipo_na_raw="${parsed_na##*|}"
                # Extrai apenas Nx se vier embutido em HxNy
                subtipo_na=$(echo "$subtipo_na_raw" | grep -oiE 'N[0-9]{1,2}' | head -1 || echo "$subtipo_na_raw")
                hit_na=$(awk 'NR==1{print $2}' "$na_blast_out")
                log_info "    NA → Subtipo: $subtipo_na | Hit: $hit_na"
            else
                log_warn "    NA → Sem hit BLAST para $sample_name."
            fi
        else
            log_warn "    NA → Arquivo vazio/ausente para $sample_name."
        fi

        echo -e "${sample_name}\t${tipo_final}\t${subtipo_ha}\t${subtipo_na}\t${hit_ha}\t${hit_na}" \
            >> "$blast_summary"
    done
    shopt -u nullglob

    log_info "Resumo BLAST: $blast_summary"
    mark_done "$checkpoint_blast"
fi

# --------------------------------------------------------------------------- #
# Etapa 5 — Nextclade: classificação de clado
# --------------------------------------------------------------------------- #
#
# v1.6.2: FluB usa dataset unificado "B" (nextstrain/flu/b/ha/KX058884).
#          O TSV gerado contém a coluna "lineage" (Victoria/Yamagata) em vez
#          de "short-clade". parse_nextclade_clade() detecta isso automaticamente.
# v1.6.1: cabeçalho do TSV de fallback inclui ambas as colunas opcionais
#          ("short-clade" e "lineage") para compatibilidade com a função.
# --------------------------------------------------------------------------- #

log_section "Etapa 5 — Nextclade: classificação de clado"

checkpoint_nextclade="$checkpoint_dir/.nextclade.done"
nextclade_dir="$output_dir/assembly_final/nextclade_results"
nextclade_datasets_dir="$DB_DIR/nextclade_datasets"
nextclade_per_sample_dir="$nextclade_dir/per_sample"
mkdir -p "$nextclade_per_sample_dir" "$nextclade_datasets_dir"

_ensure_nextclade_dataset() {
    local key="$1"
    local dataset="${NEXTCLADE_DATASETS[$key]}"
    local dest="$nextclade_datasets_dir/$key"
    if [ ! -d "$dest" ] || [ ! -f "$dest/pathogen.json" ]; then
        log_info "  Baixando dataset Nextclade '$key' ($dataset)..."
        $CONDA_RUN nextclade dataset get \
            --name "$dataset" \
            --output-dir "$dest" \
            >> "$run_log_file" 2>&1
        log_info "  Dataset '$key' pronto."
    else
        log_info "  Dataset '$key' já disponível localmente."
    fi
}

if is_done "$checkpoint_nextclade"; then
    log_skip "Nextclade já concluído (checkpoint). Pulando etapa 5."
else
    # Identifica quais datasets serão necessários antes de baixar
    datasets_needed=()
    while IFS=$'\t' read -r sample tipo subtipo_ha subtipo_na hit_ha hit_na; do
        [ "$sample" = "amostra" ] && continue
        key=$(get_nextclade_key "$tipo" "$subtipo_ha")
        if [ -n "$key" ]; then
            if ! printf '%s\n' "${datasets_needed[@]:-}" | grep -qx "$key"; then
                datasets_needed+=("$key")
            fi
        fi
    done < "$blast_summary"

    # Garante download dos datasets necessários
    for key in "${datasets_needed[@]:-}"; do
        [ -z "$key" ] && continue
        _ensure_nextclade_dataset "$key"
    done

    # Processa cada amostra
    while IFS=$'\t' read -r sample tipo subtipo_ha subtipo_na hit_ha hit_na; do
        [ "$sample" = "amostra" ] && continue

        nc_out="$nextclade_per_sample_dir/${sample}_nextclade.tsv"

        if [ -f "$nc_out" ] && [ -s "$nc_out" ]; then
            log_skip "  Nextclade já concluído para $sample."
            continue
        fi

        nc_key=$(get_nextclade_key "$tipo" "$subtipo_ha")
        ha_file="$seg_base/single_segment4/${sample}_segmento_4.fasta"

        # ── Cabeçalho de fallback ─────────────────────────────────────────
        # Inclui "short-clade" (FluA) e "lineage" (FluB) para que
        # parse_nextclade_clade() funcione corretamente em todos os casos.
        # Os valores N/A indicam ausência de dataset ou de arquivo HA.
        echo -e "seqName\tclade\tshort-clade\tlineage\tqc.overallStatus" > "$nc_out"

        if [ -z "$nc_key" ]; then
            log_warn "  $sample → Tipo '$tipo $subtipo_ha' sem dataset Nextclade (subtipo atípico). Apenas BLAST."
            echo -e "${sample}\tN/A\tN/A\tN/A\tN/A" >> "$nc_out"
            continue
        fi

        if [ ! -s "$ha_file" ]; then
            log_warn "  $sample → Arquivo HA ausente/vazio. Pulando Nextclade."
            echo -e "${sample}\tN/A\tN/A\tN/A\tN/A" >> "$nc_out"
            continue
        fi

        log_info "  Nextclade: $sample → dataset '$nc_key'"
        if ! $CONDA_RUN nextclade run \
            --input-dataset "$nextclade_datasets_dir/$nc_key" \
            --output-tsv "$nc_out" \
            "$ha_file" \
            >> "$run_log_file" 2>&1; then
            log_warn "  Nextclade falhou para $sample. Verifique $run_log_file."
            echo -e "${sample}\tERRO\tERRO\tERRO\tERRO" >> "$nc_out"
        fi

    done < "$blast_summary"

    mark_done "$checkpoint_nextclade"
    log_info "Nextclade concluído."
fi

# --------------------------------------------------------------------------- #
# Etapa 6 — Consolidação: typing_results.tsv
# --------------------------------------------------------------------------- #
#
# v1.6.2: • subtipo_HA contém apenas o componente H (ex: H5)
#         • subtipo_NA contém apenas o componente N (ex: N8)
#         • Nova coluna classificacao_BLAST (subtipo_HA + subtipo_NA → H5N8)
#         • Nova coluna segmentos (ex: 1|2|3|4|5|6|7|8 ou 2|4|8)
#           gerada a partir dos headers do multifasta em assembly_final/
# v1.6.1: leitura de colunas por nome (não por posição); fallback "—" para
#   subtipos sem dataset Nextclade (H5N1, H7N9 etc.)
# --------------------------------------------------------------------------- #

log_section "Etapa 6 — Consolidação dos resultados de tipagem"

checkpoint_consolidate="$checkpoint_dir/.consolidate.done"
final_tsv="$output_dir/assembly_final/typing_results.tsv"

if is_done "$checkpoint_consolidate"; then
    log_skip "Consolidação já realizada (checkpoint). Pulando etapa 6."
else
    echo -e "amostra\ttipo\tsubtipo_HA\tsubtipo_NA\tclassificacao_BLAST\tclado_nextclade\tqc_nextclade\tfonte_classificacao\thit_blast_HA\thit_blast_NA\tsegmentos" \
        > "$final_tsv"

    while IFS=$'\t' read -r sample tipo subtipo_ha subtipo_na hit_ha hit_na; do
        [ "$sample" = "amostra" ] && continue

        # Valores padrão: subtipo atípico sem dataset Nextclade
        clado="—"; qc="—"; fonte="BLAST"
        nc_out="$nextclade_per_sample_dir/${sample}_nextclade.tsv"

        if [ -f "$nc_out" ] && [ -s "$nc_out" ]; then

            # ── Clado: lê "clade" e "short-clade" pelo nome da coluna ────
            # Resultado: "clade/short-clade"  ex: J.2.3/2a.3a.1
            # Retorna "—" se o subtipo não tem dataset Nextclade ou se
            # o Nextclade não produziu classificação válida.
            nc_clado=$(parse_nextclade_clade "$nc_out")

            # ── QC: lê "qc.overallStatus" pelo nome da coluna ────────────
            nc_qc=$(awk -F'\t' '
                NR == 1 {
                    for (i = 1; i <= NF; i++)
                        if ($i == "qc.overallStatus") col = i
                    next
                }
                NR == 2 {
                    val = (col ? $col : "")
                    print (val == "" || val == "N/A" || val == "ERRO") ? "\xe2\x80\x94" : val
                }
            ' "$nc_out" 2>/dev/null || echo "—")

            clado="${nc_clado:-—}"
            qc="${nc_qc:-—}"

            # fonte = BLAST+Nextclade apenas quando há clado válido
            if [ "$clado" != "—" ]; then
                fonte="BLAST+Nextclade"
            fi
        fi

        # ── classificacao_BLAST: concatena subtipo_HA + subtipo_NA ────────
        # Ex: H5 + N8 → H5N8 | nd + nd → nd
        if [ "$subtipo_ha" != "nd" ] && [ "$subtipo_na" != "nd" ] \
           && [ -n "$subtipo_ha" ] && [ -n "$subtipo_na" ]; then
            classificacao_blast="${subtipo_ha}${subtipo_na}"
        elif [ "$subtipo_ha" != "nd" ] && [ -n "$subtipo_ha" ]; then
            classificacao_blast="${subtipo_ha}"
        elif [ "$subtipo_na" != "nd" ] && [ -n "$subtipo_na" ]; then
            classificacao_blast="${subtipo_na}"
        else
            classificacao_blast="nd"
        fi
        # Para Flu B (Victoria/Yamagata), usa apenas o valor de subtipo_ha
        if [ "$tipo" = "B" ]; then
            classificacao_blast="${subtipo_ha}"
        fi

        # ── Segmentos: lê headers do multifasta em assembly_final ─────────
        fasta_file="$output_dir/assembly_final/${sample}.fasta"
        if [ -f "$fasta_file" ] && [ -s "$fasta_file" ]; then
            segmentos=$(grep '^>' "$fasta_file" \
                | grep -oiE "${sample}_[0-9]+" \
                | grep -oE '[0-9]+$' \
                | sort -n \
                | tr '\n' '|' \
                | sed 's/|$//')
            [ -z "$segmentos" ] && segmentos="—"
        else
            segmentos="—"
        fi

        echo -e "${sample}\t${tipo}\t${subtipo_ha}\t${subtipo_na}\t${classificacao_blast}\t${clado}\t${qc}\t${fonte}\t${hit_ha}\t${hit_na}\t${segmentos}" \
            >> "$final_tsv"

    done < "$blast_summary"

    log_info "Tipagem consolidada: $final_tsv"
    mark_done "$checkpoint_consolidate"
fi

# --------------------------------------------------------------------------- #
# Etapa 7 — Normalização de bases em todos os FASTA
# --------------------------------------------------------------------------- #

log_section "Etapa 7 — Normalização de bases em todos os FASTA"

checkpoint_norm="$checkpoint_dir/.normalization.done"

if is_done "$checkpoint_norm"; then
    log_skip "Normalização já concluída (checkpoint). Pulando etapa 7."
else
    find "$output_dir/assembly_final/" -type f -name "*.fasta" -size +0c \
        -exec sed -i '/^>/!s/[^ACTGactg]/N/g' {} +
    log_info "Bases inválidas → 'N' em todos os FASTA."
    mark_done "$checkpoint_norm"
fi

# --------------------------------------------------------------------------- #
# Resumo final
# --------------------------------------------------------------------------- #

log_section "Pipeline finalizado"

echo ""
log_info "Resultados em: $output_dir"
log_info "  ├── assembly_final/"
log_info "  │   ├── *.fasta                  → Genomas montados"
log_info "  │   ├── typing_results.tsv        → ★ Tipagem consolidada (BLAST + Nextclade)"
log_info "  │   ├── segments/                 → Segmentos 1–8 (individuais e combinados)"
log_info "  │   ├── blast_results/            → TSVs brutos BLASTN (HA e NA)"
log_info "  │   └── nextclade_results/        → TSVs brutos Nextclade por amostra"
log_info "  └── run_log.txt                   → Log completo"
echo ""
log_info "Banco NCBI Influenza : $DB_DIR"
log_info "Modo de entrada      : $SEQ_MODE"
echo ""
log_info "Amostras OK         : $amostras_ok"
[ "$amostras_puladas" -gt 0 ] && log_info "  (retomadas         : $amostras_puladas)"
[ "$amostras_falha"   -gt 0 ] && log_warn "Amostras com falha  : $amostras_falha"
echo ""
echo "████████████████████████████████████████████████████████████████████████"
echo "█       Pipeline MK Flu-Pipe concluído — $(timestamp)                  █ "
echo "█   Desenvolvido por Jean Phellipe M. Nascimento — LACEN/AL            █"
echo "████████████████████████████████████████████████████████████████████████"
echo ""
