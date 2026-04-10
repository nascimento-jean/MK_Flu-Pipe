# MK Flu-Pipe

<div align="center">

![Version](https://img.shields.io/badge/version-1.6.2-blue)
![License](https://img.shields.io/badge/license-Academic%20%2F%20Public%20Health-green)
![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20WSL2-lightgrey)
![Language](https://img.shields.io/badge/language-Python%20%7C%20Bash-yellow)

**Ferramenta de montagem, subtipagem e análise de genomas de Influenza a partir de dados NGS**

*Genome assembly, subtyping and analysis tool for Influenza from NGS data*

Desenvolvido por / Developed by **Jean Phellipe Marques do Nascimento**
Laboratório de Vigilância Genômica — **LACEN/AL**

</div>

---

> 🇧🇷 **Português** | [🇺🇸 English below ↓](#-english-version)

---

## 🇧🇷 Versão em Português

### Índice

- [O que é o MK Flu-Pipe?](#o-que-é-o-mk-flu-pipe)
- [O que o pipeline faz?](#-o-que-o-pipeline-faz)
- [Novidades da versão 1.6.2](#-novidades-da-versão-162)
- [Pré-requisitos](#-pré-requisitos)
- [Instalação passo a passo](#️-instalação-passo-a-passo)
- [Como usar](#️-como-usar)
- [Módulos IRMA disponíveis](#-módulos-irma-disponíveis)
- [Modos de entrada (formatos de sequenciamento)](#-modos-de-entrada-formatos-de-sequenciamento)
- [Estrutura de saída](#-estrutura-de-saída)
- [Perguntas frequentes](#-perguntas-frequentes)

---

### O que é o MK Flu-Pipe?

**MK Flu-Pipe** é uma ferramenta bioinformática desenvolvida no **Laboratório de Vigilância Genômica do LACEN/AL** para montagem, subtipagem e análise de genomas do vírus Influenza a partir de dados de sequenciamento de nova geração (NGS).

A ferramenta integra três componentes principais:

- **[IRMA (CDC)](https://wonder.cdc.gov/amd/flu/irma/index.html)** — Montagem dos genomas a partir dos arquivos brutos FASTQ
- **[BLAST (NCBI)](https://blast.ncbi.nlm.nih.gov/)** — Subtipagem dos segmentos HA e NA contra o banco de dados de influenza do NCBI (atualizado automaticamente a cada 90 dias)
- **[Nextclade](https://clades.nextstrain.org/)** — Classificação de clados para H1N1pdm, H3N2, Influenza B (Victoria + Yamagata em dataset unificado) e H5Nx (todos os clados)

A ferramenta oferece tanto uma **interface gráfica (GUI)** amigável — ideal para quem não tem experiência em linha de comando — quanto uso via **linha de comando**.

---

### 🆕 Novidades da versão 1.6.2

- **Subtipo HA** agora extrai apenas o componente H (ex: `H5`, `H1`, `H3`) em vez do subtipo completo `HxNy`
- **Subtipo NA** agora extrai apenas o componente N (ex: `N1`, `N2`, `N8`)
- **Nova coluna `classificacao_BLAST`** no `typing_results.tsv`: concatena `subtipo_HA + subtipo_NA` → ex: `H5N8`, `H1N1`, `H3N2`. Para Influenza B exibe a linhagem (`Victoria` ou `Yamagata`)
- **Nova coluna `segmentos`** no `typing_results.tsv`: lista os segmentos presentes no genoma montado, separados por `|` → ex: `1|2|3|4|5|6|7|8` ou `2|4|8`
- **Suporte a H5Nx no Nextclade** via dataset comunitário (`community/moncla-lab/iav-h5/ha/all-clades`), cobrindo todos os clados H5 (H5N1, H5N8, H5N6…); exibe apenas o campo `clade` → ex: `2.3.4.4b`
- **Influenza B** agora usa dataset unificado (`nextstrain/flu/b/ha/KX058884`) que classifica Victoria e Yamagata juntos; o campo `clado_nextclade` exibe `clade/lineage` → ex: `V1A.3a.2/Victoria` ou `Y3/Yamagata`

---

### 🔬 O que o pipeline faz?

O pipeline executa as seguintes etapas automaticamente:

```
Arquivos FASTQ (R1 + R2 ou single-end)
              │
              ▼
   ┌─ Etapa 1 ─────────────────────────────────────────┐
   │  [IRMA] Montagem dos genomas por amostra          │
   └───────────────────────────────────────────────────┘
              │
              ▼
   ┌─ Etapa 2 ─────────────────────────────────────────┐
   │  Extração dos 8 segmentos genômicos               │
   │  + normalização de bases (inválidas → N)          │
   └───────────────────────────────────────────────────┘
              │
              ▼
   ┌─ Etapa 3 ─────────────────────────────────────────┐
   │  Download/atualização do banco NCBI Influenza     │
   │  (automático, a cada 90 dias)                     │
   └───────────────────────────────────────────────────┘
              │
              ▼
   ┌─ Etapa 4 ─────────────────────────────────────────┐
   │  [BLAST] Subtipagem dos segmentos HA (seg.4)      │
   │  e NA (seg.6) — determina tipo (A/B/C/D),         │
   │  componente H e componente N separadamente        │
   └───────────────────────────────────────────────────┘
              │
              ▼
   ┌─ Etapa 5 ─────────────────────────────────────────┐
   │  [Nextclade] Classificação de clado para          │
   │  H1N1pdm, H3N2, Influenza B e H5Nx               │
   │  (downloads de datasets automáticos)              │
   └───────────────────────────────────────────────────┘
              │
              ▼
   ┌─ Etapa 6 ─────────────────────────────────────────┐
   │  Consolidação de todos os resultados em           │
   │  typing_results.tsv (11 colunas)                  │
   └───────────────────────────────────────────────────┘
              │
              ▼
   ┌─ Etapa 7 ─────────────────────────────────────────┐
   │  Normalização de bases em todos os arquivos FASTA │
   │  (bases inválidas → N)                            │
   └───────────────────────────────────────────────────┘
              │
              ▼
   Resultados finais: typing_results.tsv + genomas .fasta
```

> 💡 **Sistema de checkpoint:** se o pipeline for interrompido (queda de energia, erro, etc.), ele retoma automaticamente de onde parou — sem reprocessar amostras já concluídas.

---

### 📋 Pré-requisitos

Antes de instalar o MK Flu-Pipe, você precisará instalar os seguintes componentes. Siga cada seção na ordem indicada.

> ⚠️ **Importante:** todos os comandos devem ser executados em um terminal Linux.

---

#### 1. Sistema operacional

- **Linux** — Ubuntu 20.04 ou superior (recomendado)
- **macOS** — compatível via terminal

---

#### 2. IRMA — Montador de genomas (CDC)

O IRMA é o componente responsável pela montagem dos genomas. Para instalar:

**Passo 1:** Acesse o site oficial e baixe a versão mais recente:
```
https://wonder.cdc.gov/amd/flu/irma/index.html
```
Escolha o arquivo "Universal" (compatível com Intel e ARM, Linux e macOS).

**Passo 2:** Descompacte o arquivo baixado:
```bash
# Substitua <versao> pelo número da versão baixada (ex: 1.1.5)
unzip IRMA_<versao>.zip -d ~/irma
```

**Passo 3:** Adicione o IRMA ao PATH (caminho de executáveis):
```bash
echo 'export PATH="$HOME/irma/flu-amd:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

**Passo 4:** Instale as dependências do IRMA:

```bash
# R (versão 4.0 ou superior) — necessário para o IRMA
sudo apt update
sudo apt install r-base -y
```
> O Perl 5 já vem instalado na maioria dos sistemas Linux. Não é necessário instalá-lo manualmente.

**Passo 5:** Teste se o IRMA foi instalado corretamente:
```bash
IRMA
# Deve exibir a mensagem de uso do IRMA sem erros
```

---

#### 3. Conda / Miniconda

O Conda é um gerenciador de ambientes que isola as ferramentas bioinformáticas. Se ainda não tiver o Conda instalado:

```bash
# Baixa o instalador do Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Executa o instalador (responda "yes" para todas as perguntas)
bash Miniconda3-latest-Linux-x86_64.sh

# Recarrega o terminal para aplicar as mudanças
source ~/.bashrc
```

Verifique se o Conda foi instalado:
```bash
conda --version
# Deve exibir algo como: conda 24.x.x
```

---

#### 4. Ambiente Conda `mk_flu` com BLAST e Nextclade

O pipeline usa um ambiente Conda chamado `mk_flu` que contém o BLAST e o Nextclade:

```bash
# Cria o ambiente com as ferramentas necessárias
conda create -n mk_flu -c conda-forge -c bioconda blast nextclade -y

# Aguarde a instalação (pode demorar alguns minutos)
# Verifique se funcionou:
conda run -n mk_flu blastn -version
conda run -n mk_flu nextclade --version
```

> ✅ Se ambos os comandos exibirem versões sem erros, o ambiente está pronto.

> ℹ️ O banco de dados de influenza do NCBI (`~/mk_flupipe_db/influenza.fna`) e os datasets do Nextclade (`~/mk_flupipe_db/nextclade_datasets/`) são **baixados automaticamente** na primeira execução do pipeline. Não é necessário baixá-los manualmente.

---

#### 5. Python 3 + PyGObject (somente para a interface gráfica)

A interface gráfica requer Python 3 e a biblioteca GTK. Na maioria dos sistemas Ubuntu/Debian, estas já estão instaladas. Caso não estejam:

```bash
sudo apt update
sudo apt install python3 python3-gi gir1.2-gtk-3.0 -y
```

Verifique se está funcionando:
```bash
python3 -c "import gi; print('PyGObject OK')"
```

---

### ⚙️ Instalação passo a passo

Com todos os pré-requisitos instalados, instale o MK Flu-Pipe:

**Passo 1:** Clone o repositório:
```bash
git clone https://github.com/<seu-usuario>/MK_Flu-Pipe.git
cd MK_Flu-Pipe
```

**Passo 2:** Dê permissão de execução ao script:
```bash
chmod +x script_influenza_gui.sh
```

**Passo 3:** Verifique a estrutura do projeto:
```
MK_Flu-Pipe/
├── gui_pipeline.py           ← Interface gráfica (Python/GTK)
├── script_influenza_gui.sh   ← Script principal de análise (Bash)
└── README.md                 ← Este arquivo
```

A instalação está completa!

---

### 🖥️ Como usar

#### Opção A — Interface Gráfica (recomendado para iniciantes)

Esta é a forma mais fácil de usar o pipeline. Abra um terminal na pasta do projeto e execute:

```bash
python3 gui_pipeline.py
```

Uma janela será aberta com os seguintes campos:

| Campo | O que fazer |
|-------|-------------|
| **Diretório de entrada** | Clique no botão e selecione a pasta onde estão seus arquivos FASTQ |
| **Diretório de saída** | Clique no botão e selecione (ou crie) a pasta onde os resultados serão salvos |
| **Módulo IRMA** | Escolha o módulo adequado para seus dados (veja tabela de módulos abaixo) |
| **Modo de entrada** | Deixe em "Automático" na maioria dos casos — o pipeline detecta o formato sozinho |

**Botões disponíveis:**

- **▶ Executar Pipeline** — Inicia a análise
- **⏹ Parar** — Interrompe a execução (o pipeline pode ser retomado depois)
- **💾 Salvar Log** — Salva o registro da execução em um arquivo `.txt` com timestamp
- **🗑 Limpar Log** — Limpa a área de log na tela

> 💡 A GUI salva automaticamente o último diretório, módulo e modo de entrada utilizados em `~/.mkflupipe_config.json`. Na próxima vez que abrir o programa, os campos já estarão preenchidos.

> 💡 Se fechar a janela enquanto o pipeline está rodando, um aviso será exibido perguntando se deseja interromper a execução.

---

#### Opção B — Linha de Comando

Para usuários com experiência em terminal, o pipeline pode ser executado diretamente:

**Uso básico:**
```bash
bash script_influenza_gui.sh /caminho/para/entrada /caminho/para/saida MODULO
```

**Exemplo real:**
```bash
bash script_influenza_gui.sh ~/dados/fastq ~/resultados/analise01 FLU
```

**Com modo de sequenciamento específico:**
```bash
bash script_influenza_gui.sh ~/dados/fastq ~/resultados/analise01 FLU illumina_paired
```

**Modo interativo** (o script pergunta os parâmetros um por um):
```bash
bash script_influenza_gui.sh
```

---

### 🦠 Módulos IRMA disponíveis

| Módulo | Quando usar |
|--------|-------------|
| `FLU` | **Uso geral — recomendado para a maioria dos casos.** Suporta Influenza A e B de dados Illumina paired-end. |
| `FLU-utr` | Quando você precisa das regiões UTR (não-codificantes) incluídas nas sequências consenso. |
| `FLU-lowQC` | Para amostras com baixa qualidade de sequenciamento (ex: amostras degradadas ou com pouco material). |
| `FLU_AD` | Quando há suspeita ou necessidade de detectar também Influenza C e D, além de A e B. |
| `FLU-minion` | Para dados gerados pelo sequenciador Oxford Nanopore (leituras longas). |

---

### 📂 Modos de entrada (formatos de sequenciamento)

O pipeline detecta automaticamente o formato dos seus arquivos FASTQ na seguinte ordem de prioridade:

| Modo | Padrão de arquivo esperado | Origem típica |
|------|---------------------------|---------------|
| `illumina_paired` | `amostra_L001_R1_001.fastq.gz` + `R2` | Illumina BaseSpace |
| `sra_paired` | `amostra_1.fastq.gz` + `amostra_2.fastq.gz` | Download NCBI SRA |
| `generic_paired` | `amostra_R1_001.fastq.gz` + `R2` | Illumina genérico |
| `single` | `amostra.fastq.gz` (um arquivo por amostra) | ONT, Ion Torrent, PacBio, SRA single |

> ℹ️ Em caso de dúvida, deixe o modo em **Automático**. O pipeline verificará os padrões de nome dos arquivos e detectará o formato correto.

---

### 📁 Estrutura de saída

Após a execução bem-sucedida, o diretório de saída conterá:

```
saida/
│
├── <nome_amostra>/                    ← Saída bruta do IRMA por amostra
│   ├── amended_consensus/             ← Sequências consenso do IRMA
│   ├── tables/                        ← Tabelas de variantes
│   └── ...
│
├── assembly_final/
│   ├── <amostra>.fasta                ← Genoma completo concatenado (todos os segmentos)
│   │
│   ├── segments/
│   │   ├── segmento_1.fasta           ← Segmento 1 de todas as amostras combinadas
│   │   ├── segmento_2.fasta
│   │   ├── ...
│   │   ├── segmento_8.fasta
│   │   ├── single_segment4/           ← HA (seg.4) — um arquivo FASTA por amostra
│   │   ├── single_segment6/           ← NA (seg.6) — um arquivo FASTA por amostra
│   │   └── single_segment7/           ← M  (seg.7) — um arquivo FASTA por amostra
│   │
│   ├── blast_results/                 ← TSVs brutos do BLASTN (HA e NA)
│   ├── nextclade_results/             ← TSVs brutos do Nextclade por amostra
│   └── typing_results.tsv             ← ★ RESULTADO PRINCIPAL: tipagem consolidada
│
├── run_log.txt                        ← Log completo de execução
└── .checkpoints/                      ← Arquivos de controle do sistema de checkpoint
```

### Entendendo o arquivo `typing_results.tsv`

Este é o arquivo mais importante da saída. Ele contém uma linha por amostra com as seguintes colunas:

| Coluna | Descrição |
|--------|-----------|
| `amostra` | Nome da amostra |
| `tipo` | Tipo de Influenza: A, B, C ou D |
| `subtipo_HA` | Componente H do segmento HA (ex: `H1`, `H3`, `H5`) |
| `subtipo_NA` | Componente N do segmento NA (ex: `N1`, `N2`, `N8`) |
| `classificacao_BLAST` | Subtipo completo concatenado (ex: `H1N1`, `H3N2`, `H5N8`). Para Flu B exibe a linhagem (`Victoria` ou `Yamagata`) |
| `clado_nextclade` | Clado / sub-clado: `clade/short-clade` para FluA (ex: `J.2.3/2a.3a.1`), `clade/lineage` para FluB (ex: `V1A.3a.2/Victoria`), apenas `clade` para H5Nx (ex: `2.3.4.4b`), ou `—` se sem dataset |
| `qc_nextclade` | Status de qualidade do Nextclade (`good`, `mediocre`, `bad`, ou `—`) |
| `fonte_classificacao` | `BLAST` ou `BLAST+Nextclade` dependendo se o clado foi determinado |
| `hit_blast_HA` | Melhor hit BLAST para o segmento HA |
| `hit_blast_NA` | Melhor hit BLAST para o segmento NA |
| `segmentos` | Segmentos presentes no genoma montado, separados por `\|` (ex: `1\|2\|3\|4\|5\|6\|7\|8` ou `2\|4\|8`) |

---

### Datasets Nextclade disponíveis

| Chave interna | Dataset | Subtipos cobertos |
|---------------|---------|-------------------|
| `H1N1` | `nextstrain/flu/h1n1pdm/ha/MW626062` | H1N1pdm |
| `H3N2` | `nextstrain/flu/h3n2/ha/EPI1857216` | H3N2 |
| `B` | `nextstrain/flu/b/ha/KX058884` | Influenza B (Victoria + Yamagata) |
| `H5` | `community/moncla-lab/iav-h5/ha/all-clades` | Todos os clados H5Nx |

> ℹ️ Os datasets são baixados automaticamente na primeira execução e armazenados em `~/mk_flupipe_db/nextclade_datasets/`. Subtipos sem dataset (ex: H7N9, H9N2) recebem `—` nas colunas de clado e `fonte_classificacao = BLAST`.

---

### ❓ Perguntas frequentes

**Os arquivos FASTQ precisam seguir algum padrão de nome?**
Sim. O pipeline precisa identificar o nome da amostra a partir do nome do arquivo. O modo `illumina_paired` espera o padrão Illumina BaseSpace: `<amostra>_L001_R1_001.fastq.gz` e `<amostra>_L001_R2_001.fastq.gz`. Para outros formatos, use o modo de detecção automática ou selecione o modo manualmente na GUI.

**E se uma amostra não gerar resultado no IRMA?**
O pipeline registra a falha no `run_log.txt`, pula a amostra e continua processando as demais. Ao final, o resumo mostrará quantas amostras foram processadas com sucesso e quantas falharam.

**O pipeline foi interrompido no meio. Preciso recomeçar do zero?**
Não. O MK Flu-Pipe possui um **sistema de checkpoint** — ao ser executado novamente com os mesmos diretórios de entrada e saída, ele detecta automaticamente quais amostras já foram processadas e retoma de onde parou.

**Posso usar no Windows?**
Sim. Instale o WSL2 (Windows Subsystem for Linux) e siga as instruções normais de instalação Linux dentro do WSL2. No Windows 11 com WSLg habilitado, a interface gráfica (GUI) também funciona diretamente, sem configuração adicional.

**O pipeline funciona com dados Nanopore (Oxford Nanopore)?**
Sim. Use o módulo `FLU-minion` e selecione o modo de entrada `single` (ou deixe em automático).

**O que significa `—` nos campos de clado no resultado?**
Significa que o Nextclade não possui um dataset disponível para aquele subtipo específico (ex: H7N9, H9N2, H2N2). Nestes casos, a tipagem é feita apenas pelo BLAST e a coluna `fonte_classificacao` mostrará apenas `BLAST`.

**O que é a coluna `classificacao_BLAST`?**
É o subtipo completo montado a partir dos resultados BLAST, combinando os componentes H e N identificados separadamente: `subtipo_HA + subtipo_NA` → ex: `H5N8`, `H1N1`. Para Influenza B, exibe a linhagem (`Victoria` ou `Yamagata`). Esta coluna facilita a leitura rápida do subtipo sem precisar combinar as colunas `subtipo_HA` e `subtipo_NA` manualmente.

**O que é a coluna `segmentos`?**
Lista os segmentos genômicos presentes no genoma montado para aquela amostra, separados por `|`. Um valor `1|2|3|4|5|6|7|8` indica genoma completo; valores como `2|4|8` indicam que apenas esses segmentos foram montados, o que pode ocorrer em amostras com baixa cobertura.

**Quanto tempo leva a execução?**
Depende do número de amostras e do hardware disponível. O pipeline usa automaticamente todos os núcleos de CPU disponíveis. Para uma corrida típica de 12–24 amostras Illumina, espere entre 30 minutos e 2 horas.

**Onde ficam os datasets do Nextclade e o banco NCBI?**
Ambos são armazenados em `~/mk_flupipe_db/`: o banco NCBI Influenza em `influenza.fna` e os datasets do Nextclade em `nextclade_datasets/`. Eles são reutilizados em execuções subsequentes. O banco NCBI é atualizado automaticamente a cada 90 dias.

---

### 👨‍🔬 Autoria e contato

Desenvolvido por **Jean Phellipe Marques do Nascimento**
Laboratório de Vigilância Genômica — **LACEN/AL** (Laboratório Central de Saúde Pública de Alagoas)

Para dúvidas, problemas ou sugestões, abra uma [**Issue**](https://github.com/nascimento-jean/MK_Flu-Pipe/issues) neste repositório.

---

### 📄 Licença

Este projeto é distribuído para uso acadêmico e de saúde pública.
Consulte as licenças das ferramentas integradas:
- [Licença do IRMA (CDC)](https://wonder.cdc.gov/amd/flu/irma/disclaimer.html)
- [Licença do BLAST (NCBI)](https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/scripts/projects/blast/LICENSE)
- [Licença do Nextclade (Nextstrain)](https://github.com/nextstrain/nextclade/blob/master/LICENSE)

---
---

<a name="-english-version"></a>
## 🇺🇸 English Version

### Table of Contents

- [What is MK Flu-Pipe?](#what-is-mk-flu-pipe)
- [What does the pipeline do?](#-what-does-the-pipeline-do)
- [What's new in v1.6.2](#-whats-new-in-v162)
- [Prerequisites](#-prerequisites)
- [Step-by-step installation](#️-step-by-step-installation)
- [How to use](#️-how-to-use)
- [Available IRMA modules](#-available-irma-modules)
- [Input modes (sequencing formats)](#-input-modes-sequencing-formats)
- [Output structure](#-output-structure)
- [FAQ](#-faq)

---

### What is MK Flu-Pipe?

**MK Flu-Pipe** is a bioinformatics tool developed at the **Genomic Surveillance Laboratory of LACEN/AL** (Alagoas Central Public Health Laboratory, Brazil) for the assembly, subtyping, and analysis of Influenza virus genomes from next-generation sequencing (NGS) data.

The tool integrates three main components:

- **[IRMA (CDC)](https://wonder.cdc.gov/amd/flu/irma/index.html)** — Genome assembly from raw FASTQ files
- **[BLAST (NCBI)](https://blast.ncbi.nlm.nih.gov/)** — Subtyping of HA and NA segments against the NCBI Influenza database (automatically updated every 90 days)
- **[Nextclade](https://clades.nextstrain.org/)** — Clade classification for H1N1pdm, H3N2, Influenza B (Victoria + Yamagata via unified dataset), and H5Nx (all clades)

The tool provides both a user-friendly **graphical interface (GUI)** — ideal for users without command-line experience — and **command-line** usage.

---

### 🆕 What's new in v1.6.2

- **HA subtype** now extracts only the H component (e.g., `H5`, `H1`, `H3`) instead of the full `HxNy` subtype
- **NA subtype** now extracts only the N component (e.g., `N1`, `N2`, `N8`)
- **New `classificacao_BLAST` column** in `typing_results.tsv`: concatenates `subtipo_HA + subtipo_NA` → e.g., `H5N8`, `H1N1`, `H3N2`. For Influenza B, displays the lineage (`Victoria` or `Yamagata`)
- **New `segmentos` column** in `typing_results.tsv`: lists the genomic segments present in the assembled genome, separated by `|` → e.g., `1|2|3|4|5|6|7|8` or `2|4|8`
- **H5Nx support in Nextclade** via community dataset (`community/moncla-lab/iav-h5/ha/all-clades`), covering all H5 clades (H5N1, H5N8, H5N6…); displays only the `clade` field → e.g., `2.3.4.4b`
- **Influenza B** now uses a unified dataset (`nextstrain/flu/b/ha/KX058884`) that classifies Victoria and Yamagata together; the `clado_nextclade` field shows `clade/lineage` → e.g., `V1A.3a.2/Victoria` or `Y3/Yamagata`

---

### 🔬 What does the pipeline do?

The pipeline runs the following steps automatically:

```
FASTQ files (R1 + R2 or single-end)
              │
              ▼
   ┌─ Step 1 ──────────────────────────────────────────┐
   │  [IRMA] Genome assembly per sample               │
   └───────────────────────────────────────────────────┘
              │
              ▼
   ┌─ Step 2 ──────────────────────────────────────────┐
   │  Extraction of all 8 genomic segments            │
   │  + base normalization (invalid bases → N)        │
   └───────────────────────────────────────────────────┘
              │
              ▼
   ┌─ Step 3 ──────────────────────────────────────────┐
   │  Download/update of the NCBI Influenza database   │
   │  (automatic, every 90 days)                       │
   └───────────────────────────────────────────────────┘
              │
              ▼
   ┌─ Step 4 ──────────────────────────────────────────┐
   │  [BLAST] Subtyping of HA (seg.4) and NA (seg.6)  │
   │  — determines type (A/B/C/D) and H/N components  │
   │  separately (e.g. H5 and N8)                     │
   └───────────────────────────────────────────────────┘
              │
              ▼
   ┌─ Step 5 ──────────────────────────────────────────┐
   │  [Nextclade] Clade classification for H1N1pdm,   │
   │  H3N2, Influenza B, and H5Nx                     │
   │  (datasets downloaded automatically)             │
   └───────────────────────────────────────────────────┘
              │
              ▼
   ┌─ Step 6 ──────────────────────────────────────────┐
   │  Consolidation of all results into               │
   │  typing_results.tsv (11 columns)                 │
   └───────────────────────────────────────────────────┘
              │
              ▼
   ┌─ Step 7 ──────────────────────────────────────────┐
   │  Base normalization across all FASTA files       │
   │  (invalid bases → N)                             │
   └───────────────────────────────────────────────────┘
              │
              ▼
   Final results: typing_results.tsv + .fasta genomes
```

> 💡 **Checkpoint system:** if the pipeline is interrupted (power failure, error, etc.), it automatically resumes from where it stopped — without reprocessing already-completed samples.

---

### 📋 Prerequisites

Before installing MK Flu-Pipe, you need to install the following components. Follow each section in order.

> ⚠️ **Important:** all commands must be run in a Linux terminal.
---

#### 1. Operating System

- **Linux** — Ubuntu 20.04 or later (recommended)
- **macOS** — compatible via terminal

---

#### 2. IRMA — Genome Assembler (CDC)

IRMA is the core assembly component. To install:

**Step 1:** Visit the official website and download the latest version:
```
https://wonder.cdc.gov/amd/flu/irma/index.html
```
Choose the "Universal" file (compatible with Intel and ARM, Linux and macOS).

**Step 2:** Unzip the downloaded file:
```bash
# Replace <version> with the downloaded version number (e.g., 1.1.5)
unzip IRMA_<version>.zip -d ~/irma
```

**Step 3:** Add IRMA to your PATH:
```bash
echo 'export PATH="$HOME/irma/flu-amd:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

**Step 4:** Install IRMA dependencies:
```bash
# R (version 4.0 or later) — required for IRMA
sudo apt update
sudo apt install r-base -y
```
> Perl 5 is already installed on most Linux systems and does not need to be installed manually.

**Step 5:** Test that IRMA was installed correctly:
```bash
IRMA
# Should display IRMA usage information without errors
```

---

#### 3. Conda / Miniconda

Conda is an environment manager that isolates bioinformatics tools. If you do not have Conda installed yet:

```bash
# Download the Miniconda installer
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Run the installer (answer "yes" to all prompts)
bash Miniconda3-latest-Linux-x86_64.sh

# Reload the terminal to apply changes
source ~/.bashrc
```

Verify the installation:
```bash
conda --version
# Should display something like: conda 24.x.x
```

---

#### 4. `mk_flu` Conda environment with BLAST and Nextclade

The pipeline uses a Conda environment named `mk_flu` containing BLAST and Nextclade:

```bash
# Create the environment with the required tools
conda create -n mk_flu -c conda-forge -c bioconda blast nextclade -y

# Wait for installation (may take a few minutes)
# Verify it worked:
conda run -n mk_flu blastn -version
conda run -n mk_flu nextclade --version
```

> ✅ If both commands display version numbers without errors, the environment is ready.

> ℹ️ The NCBI Influenza database (`~/mk_flupipe_db/influenza.fna`) and Nextclade datasets (`~/mk_flupipe_db/nextclade_datasets/`) are **downloaded automatically** during the first pipeline run. No manual download is required.

---

#### 5. Python 3 + PyGObject (GUI only)

The graphical interface requires Python 3 and the GTK library. On most Ubuntu/Debian systems, these are already installed. If not:

```bash
sudo apt update
sudo apt install python3 python3-gi gir1.2-gtk-3.0 -y
```

Verify it is working:
```bash
python3 -c "import gi; print('PyGObject OK')"
```

---

### ⚙️ Step-by-step installation

With all prerequisites installed, install MK Flu-Pipe:

**Step 1:** Clone the repository:
```bash
git clone https://github.com/<your-username>/MK_Flu-Pipe.git
cd MK_Flu-Pipe
```

**Step 2:** Make the script executable:
```bash
chmod +x script_influenza_gui.sh
```

**Step 3:** Verify the project structure:
```
MK_Flu-Pipe/
├── gui_pipeline.py           ← Graphical interface (Python/GTK)
├── script_influenza_gui.sh   ← Main analysis script (Bash)
└── README.md                 ← This file
```

Installation is complete!

---

### 🖥️ How to use

#### Option A — Graphical Interface (recommended for beginners)

This is the easiest way to run the pipeline. Open a terminal in the project folder and run:

```bash
python3 gui_pipeline.py
```

A window will open with the following fields:

| Field | What to do |
|-------|------------|
| **Input directory** | Click the button and select the folder containing your FASTQ files |
| **Output directory** | Click the button and select (or create) the folder where results will be saved |
| **IRMA module** | Choose the appropriate module for your data (see module table below) |
| **Input mode** | Leave as "Automatic" in most cases — the pipeline detects the format on its own |

**Available buttons:**

- **▶ Run Pipeline** — Starts the analysis
- **⏹ Stop** — Interrupts the run (the pipeline can be resumed later)
- **💾 Save Log** — Saves the execution log to a `.txt` file with a timestamp
- **🗑 Clear Log** — Clears the log area on screen

> 💡 The GUI automatically saves the last directory, module, and input mode used in `~/.mkflupipe_config.json`. Next time you open the program, the fields will already be filled in.

> 💡 If you close the window while the pipeline is running, a warning will ask whether you want to interrupt the execution.

---

#### Option B — Command Line

For users experienced with the terminal, the pipeline can be run directly:

**Basic usage:**
```bash
bash script_influenza_gui.sh /path/to/input /path/to/output MODULE
```

**Real example:**
```bash
bash script_influenza_gui.sh ~/data/fastq ~/results/analysis01 FLU
```

**With a specific sequencing mode:**
```bash
bash script_influenza_gui.sh ~/data/fastq ~/results/analysis01 FLU illumina_paired
```

**Interactive mode** (the script asks for parameters one by one):
```bash
bash script_influenza_gui.sh
```

---

### 🦠 Available IRMA modules

| Module | When to use |
|--------|-------------|
| `FLU` | **General use — recommended for most cases.** Supports Influenza A and B from Illumina paired-end data. |
| `FLU-utr` | When you need UTR (non-coding) regions included in the consensus sequences. |
| `FLU-lowQC` | For samples with low sequencing quality (e.g., degraded samples or low input material). |
| `FLU_AD` | When there is a need to detect Influenza C and D in addition to A and B. |
| `FLU-minion` | For data generated by Oxford Nanopore sequencers (long reads). |

---

### 📂 Input modes (sequencing formats)

The pipeline automatically detects your FASTQ file format in the following priority order:

| Mode | Expected file pattern | Typical origin |
|------|-----------------------|----------------|
| `illumina_paired` | `sample_L001_R1_001.fastq.gz` + `R2` | Illumina BaseSpace |
| `sra_paired` | `sample_1.fastq.gz` + `sample_2.fastq.gz` | NCBI SRA download |
| `generic_paired` | `sample_R1_001.fastq.gz` + `R2` | Generic Illumina |
| `single` | `sample.fastq.gz` (one file per sample) | ONT, Ion Torrent, PacBio, SRA single |

> ℹ️ When in doubt, leave the mode as **Automatic**. The pipeline checks file naming patterns and detects the correct format.

---

### 📁 Output structure

After a successful run, the output directory will contain:

```
output/
│
├── <sample_name>/                     ← Raw IRMA output per sample
│   ├── amended_consensus/             ← IRMA consensus sequences
│   ├── tables/                        ← Variant tables
│   └── ...
│
├── assembly_final/
│   ├── <sample>.fasta                 ← Full concatenated genome (all segments)
│   │
│   ├── segments/
│   │   ├── segmento_1.fasta           ← Segment 1 from all samples combined
│   │   ├── segmento_2.fasta
│   │   ├── ...
│   │   ├── segmento_8.fasta
│   │   ├── single_segment4/           ← HA (seg.4) — one FASTA per sample
│   │   ├── single_segment6/           ← NA (seg.6) — one FASTA per sample
│   │   └── single_segment7/           ← M  (seg.7) — one FASTA per sample
│   │
│   ├── blast_results/                 ← Raw BLASTN TSVs (HA and NA)
│   ├── nextclade_results/             ← Raw Nextclade TSVs per sample
│   └── typing_results.tsv             ← ★ MAIN RESULT: consolidated typing
│
├── run_log.txt                        ← Full execution log
└── .checkpoints/                      ← Checkpoint control files
```

### Understanding `typing_results.tsv`

This is the most important output file. It contains one row per sample with the following columns:

| Column | Description |
|--------|-------------|
| `amostra` | Sample name |
| `tipo` | Influenza type: A, B, C, or D |
| `subtipo_HA` | H component of the HA segment (e.g., `H1`, `H3`, `H5`) |
| `subtipo_NA` | N component of the NA segment (e.g., `N1`, `N2`, `N8`) |
| `classificacao_BLAST` | Full subtype assembled from BLAST (e.g., `H1N1`, `H3N2`, `H5N8`). For Flu B, displays the lineage (`Victoria` or `Yamagata`) |
| `clado_nextclade` | Clade/sub-clade: `clade/short-clade` for FluA (e.g., `J.2.3/2a.3a.1`), `clade/lineage` for FluB (e.g., `V1A.3a.2/Victoria`), `clade` only for H5Nx (e.g., `2.3.4.4b`), or `—` if no dataset |
| `qc_nextclade` | Nextclade quality status (`good`, `mediocre`, `bad`, or `—`) |
| `fonte_classificacao` | `BLAST` or `BLAST+Nextclade` depending on whether a clade was determined |
| `hit_blast_HA` | Best BLAST hit for the HA segment |
| `hit_blast_NA` | Best BLAST hit for the NA segment |
| `segmentos` | Genomic segments present in the assembled genome, separated by `\|` (e.g., `1\|2\|3\|4\|5\|6\|7\|8` or `2\|4\|8`) |

---

### Available Nextclade Datasets

| Internal key | Dataset | Subtypes covered |
|--------------|---------|-----------------|
| `H1N1` | `nextstrain/flu/h1n1pdm/ha/MW626062` | H1N1pdm |
| `H3N2` | `nextstrain/flu/h3n2/ha/EPI1857216` | H3N2 |
| `B` | `nextstrain/flu/b/ha/KX058884` | Influenza B (Victoria + Yamagata) |
| `H5` | `community/moncla-lab/iav-h5/ha/all-clades` | All H5Nx clades |

> ℹ️ Datasets are downloaded automatically on the first run and stored in `~/mk_flupipe_db/nextclade_datasets/`. Subtypes without a dataset (e.g., H7N9, H9N2) receive `—` in the clade columns and `fonte_classificacao = BLAST`.

---

### ❓ FAQ

**Do FASTQ files need to follow a specific naming convention?**
Yes. The pipeline identifies sample names from file names. The `illumina_paired` mode expects the standard Illumina BaseSpace pattern: `<sample>_L001_R1_001.fastq.gz` and `<sample>_L001_R2_001.fastq.gz`. For other formats, use automatic detection or manually select the mode in the GUI.

**What if a sample fails in IRMA?**
The pipeline logs the failure in `run_log.txt`, skips the sample, and continues processing the remaining ones. At the end, the summary shows how many samples were successfully processed and how many failed.

**The pipeline was interrupted halfway. Do I need to start over?**
No. MK Flu-Pipe has a **checkpoint system** — when rerun with the same input and output directories, it automatically detects which samples were already processed and resumes from where it stopped.

**Can I use this on Windows?**
Yes. Install WSL2 (Windows Subsystem for Linux) and follow the standard Linux installation instructions inside WSL2. On Windows 11 with WSLg enabled, the graphical interface (GUI) also works natively without additional configuration.

**Does the pipeline work with Nanopore (Oxford Nanopore) data?**
Yes. Use the `FLU-minion` module and select the `single` input mode (or leave it on automatic).

**What does `—` mean in the clade fields of the results?**
It means Nextclade does not have an available dataset for that specific subtype (e.g., H7N9, H9N2, H2N2). In these cases, typing is performed by BLAST only and the `fonte_classificacao` column will show just `BLAST`.

**What is the `classificacao_BLAST` column?**
It is the full subtype assembled from BLAST results, combining the H and N components identified separately: `subtipo_HA + subtipo_NA` → e.g., `H5N8`, `H1N1`. For Influenza B, it displays the lineage (`Victoria` or `Yamagata`). This column makes it easy to read the subtype at a glance without manually combining `subtipo_HA` and `subtipo_NA`.

**What is the `segmentos` column?**
It lists the genomic segments present in the assembled genome for that sample, separated by `|`. A value of `1|2|3|4|5|6|7|8` indicates a complete genome; values like `2|4|8` indicate that only those segments were assembled, which can occur in samples with low sequencing coverage.

**How long does a run take?**
It depends on the number of samples and available hardware. The pipeline automatically uses all available CPU cores. For a typical run of 12–24 Illumina samples, expect between 30 minutes and 2 hours.

**Where are the Nextclade datasets and NCBI database stored?**
Both are stored in `~/mk_flupipe_db/`: the NCBI Influenza database as `influenza.fna` and Nextclade datasets in `nextclade_datasets/`. They are reused in subsequent runs. The NCBI database is automatically updated every 90 days.

---

### 👨‍🔬 Authorship and contact

Developed by **Jean Phellipe Marques do Nascimento**
Genomic Surveillance Laboratory — **LACEN/AL** (Alagoas Central Public Health Laboratory, Brazil)

For questions, issues, or suggestions, please open an [**Issue**](https://github.com/nascimento-jean/MK_Flu-Pipe/issues) in this repository.

---

### 📄 License

This project is distributed for academic and public health use.
Please refer to the licenses of the integrated tools:
- [IRMA License (CDC)](https://wonder.cdc.gov/amd/flu/irma/disclaimer.html)
- [BLAST License (NCBI)](https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/scripts/projects/blast/LICENSE)
- [Nextclade License (Nextstrain)](https://github.com/nextstrain/nextclade/blob/master/LICENSE)
