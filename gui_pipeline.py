#!/usr/bin/env python3
# =============================================================================
# MK Flu-Pipe — Interface Gráfica
# Versão: 1.6.2
# Desenvolvido por: Jean Phellipe Marques do Nascimento
# Laboratório de Vigilância Genômica — LACEN/AL
# =============================================================================

import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, GLib, Pango

import threading
import subprocess
import os
import datetime
import json

# Arquivo para persistência das últimas configurações usadas
CONFIG_FILE = os.path.join(os.path.expanduser("~"), ".mkflupipe_config.json")


def load_config() -> dict:
    """Carrega as últimas configurações salvas pelo usuário."""
    try:
        with open(CONFIG_FILE, "r") as f:
            return json.load(f)
    except Exception:
        return {}


def save_config(data: dict):
    """Salva as configurações atuais para a próxima sessão."""
    try:
        with open(CONFIG_FILE, "w") as f:
            json.dump(data, f, indent=2)
    except Exception:
        pass


class JanelaPrincipal(Gtk.Window):
    def __init__(self):
        Gtk.Window.__init__(self, title="MK Flu‑Pipe by Jean Nascimento v1.6.2")
        self.set_border_width(12)
        self.set_default_size(950, 680)
        self.set_resizable(True)
        self.connect("delete-event", self._on_window_delete)

        self._processo = None       # referência ao subprocess ativo
        self._executando = False    # flag de controle de execução
        self._config = load_config()

        # ── Layout principal ───────────────────────────────────────────────
        vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=10)
        vbox.set_margin_top(8)
        vbox.set_margin_bottom(8)
        vbox.set_margin_start(8)
        vbox.set_margin_end(8)
        self.add(vbox)

        # ── Cabeçalho ─────────────────────────────────────────────────────
        header_label = Gtk.Label()
        header_label.set_markup(
            "<span size='large' weight='bold'>MK Flu‑Pipe</span>  "
            "<span size='small' foreground='gray'>v1.6.2</span>"
        )
        header_label.set_halign(Gtk.Align.START)
        vbox.pack_start(header_label, False, False, 0)

        separator = Gtk.Separator(orientation=Gtk.Orientation.HORIZONTAL)
        vbox.pack_start(separator, False, False, 2)

        # ── Grade de campos ────────────────────────────────────────────────
        grid = Gtk.Grid(column_spacing=12, row_spacing=10)
        grid.set_margin_top(4)
        vbox.pack_start(grid, False, False, 0)

        # Campo: Diretório de entrada
        lbl_in = Gtk.Label(label="Diretório de entrada:")
        lbl_in.set_halign(Gtk.Align.END)
        self.fc_input = Gtk.FileChooserButton(
            title="Selecionar diretório de entrada",
            action=Gtk.FileChooserAction.SELECT_FOLDER)
        self.fc_input.set_tooltip_text(
            "Pasta contendo os arquivos FASTQ (.fastq.gz)\n"
            "O modo de entrada é detectado automaticamente conforme o padrão dos arquivos.\n"
            "Você também pode forçar um modo específico no campo abaixo."
        )
        self.fc_input.set_hexpand(True)
        # Restaura último diretório de entrada
        if self._config.get("input_dir") and os.path.isdir(self._config["input_dir"]):
            self.fc_input.set_filename(self._config["input_dir"])

        # Campo: Diretório de saída
        lbl_out = Gtk.Label(label="Diretório de saída:")
        lbl_out.set_halign(Gtk.Align.END)
        self.fc_output = Gtk.FileChooserButton(
            title="Selecionar diretório de saída",
            action=Gtk.FileChooserAction.SELECT_FOLDER)
        self.fc_output.set_tooltip_text(
            "Pasta onde os resultados serão salvos.\n"
            "Se não existir, será criada automaticamente."
        )
        self.fc_output.set_hexpand(True)
        if self._config.get("output_dir") and os.path.isdir(self._config["output_dir"]):
            self.fc_output.set_filename(self._config["output_dir"])

        # Campo: Módulo IRMA
        lbl_mod = Gtk.Label(label="Módulo IRMA:")
        lbl_mod.set_halign(Gtk.Align.END)
        self.combo_module = Gtk.ComboBoxText()
        modules = ["FLU", "FLU-utr", "FLU-lowQC", "FLU_AD", "FLU-minion"]
        for mod in modules:
            self.combo_module.append_text(mod)
        # Restaura último módulo selecionado
        last_module = self._config.get("irma_module", "FLU")
        idx = modules.index(last_module) if last_module in modules else 0
        self.combo_module.set_active(idx)
        self.combo_module.set_tooltip_text(
            "FLU       → Influenza A e B (uso geral)\n"
            "FLU-utr   → Inclui regiões UTR\n"
            "FLU-lowQC → Para amostras de baixa qualidade\n"
            "FLU_AD    → Influenza A, B, C e D\n"
            "FLU-minion → Dados de sequenciamento Oxford Nanopore"
        )

        # Campo: Modo de sequenciamento
        lbl_mode = Gtk.Label(label="Modo de entrada:")
        lbl_mode.set_halign(Gtk.Align.END)
        self.combo_mode = Gtk.ComboBoxText()
        modes = [
            ("auto",            "Automático (detectar pelo padrão dos arquivos)"),
            ("illumina_paired", "Illumina paired-end  (*_L001_R1_001.fastq.gz)"),
            ("sra_paired",      "SRA/NCBI paired-end  (*_1.fastq.gz + *_2.fastq.gz)"),
            ("generic_paired",  "Paired-end genérico  (*_R1_*.fastq.gz)"),
            ("single",          "Single-end  (ONT / Ion Torrent / PacBio / SRA single)"),
        ]
        for mode_id, mode_label in modes:
            self.combo_mode.append(mode_id, mode_label)
        last_mode = self._config.get("seq_mode", "auto")
        self.combo_mode.set_active_id(last_mode if last_mode in dict(modes) else "auto")
        self.combo_mode.set_tooltip_text(
            "Automático     → O script detecta o formato pelos nomes dos arquivos\n"
            "Illumina       → *_L001_R1_001.fastq.gz / *_L001_R2_001.fastq.gz\n"
            "SRA paired     → *_1.fastq.gz + *_2.fastq.gz (download do NCBI SRA)\n"
            "Paired genérico→ *_R1_*.fastq.gz + *_R2_*.fastq.gz\n"
            "Single-end     → um único .fastq.gz por amostra (ONT, Ion Torrent, SRA single)"
        )

        # Posiciona na grade
        grid.attach(lbl_in,            0, 0, 1, 1)
        grid.attach(self.fc_input,     1, 0, 3, 1)

        grid.attach(lbl_out,           0, 1, 1, 1)
        grid.attach(self.fc_output,    1, 1, 3, 1)

        grid.attach(lbl_mod,           0, 2, 1, 1)
        grid.attach(self.combo_module, 1, 2, 1, 1)

        grid.attach(lbl_mode,          0, 3, 1, 1)
        grid.attach(self.combo_mode,   1, 3, 3, 1)

        # ── Botões de ação ─────────────────────────────────────────────────
        btn_box = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
        btn_box.set_margin_top(4)

        self.btn_run = Gtk.Button(label="▶  Executar Pipeline")
        self.btn_run.connect("clicked", self.on_run_clicked)
        self.btn_run.get_style_context().add_class("suggested-action")
        self.btn_run.set_tooltip_text("Inicia a montagem e tipagem das amostras")

        self.btn_stop = Gtk.Button(label="⏹  Parar")
        self.btn_stop.connect("clicked", self.on_stop_clicked)
        self.btn_stop.set_sensitive(False)
        self.btn_stop.set_tooltip_text("Interrompe a execução do pipeline")

        self.btn_save_log = Gtk.Button(label="💾  Salvar Log")
        self.btn_save_log.connect("clicked", self.on_save_log_clicked)
        self.btn_save_log.set_tooltip_text("Salva o conteúdo do log em um arquivo de texto")

        self.btn_clear_log = Gtk.Button(label="🗑  Limpar Log")
        self.btn_clear_log.connect("clicked", self.on_clear_log_clicked)
        self.btn_clear_log.set_tooltip_text("Limpa a área de log")

        btn_box.pack_start(self.btn_run,      False, False, 0)
        btn_box.pack_start(self.btn_stop,     False, False, 0)
        btn_box.pack_end(self.btn_clear_log,  False, False, 0)
        btn_box.pack_end(self.btn_save_log,   False, False, 0)

        vbox.pack_start(btn_box, False, False, 0)

        separator2 = Gtk.Separator(orientation=Gtk.Orientation.HORIZONTAL)
        vbox.pack_start(separator2, False, False, 2)

        # ── Área de log ────────────────────────────────────────────────────
        log_label = Gtk.Label(label="Log de execução:")
        log_label.set_halign(Gtk.Align.START)
        vbox.pack_start(log_label, False, False, 0)

        self.scrolled = Gtk.ScrolledWindow()
        self.scrolled.set_policy(Gtk.PolicyType.AUTOMATIC, Gtk.PolicyType.AUTOMATIC)
        self.textview = Gtk.TextView()
        self.textview.set_editable(False)
        self.textview.set_wrap_mode(Gtk.WrapMode.WORD_CHAR)
        self.textview.set_monospace(True)
        # Fonte monoespaçada para melhor legibilidade de logs
        font_desc = Pango.FontDescription("Monospace 9")
        self.textview.override_font(font_desc)
        self.buffer = self.textview.get_buffer()
        self.scrolled.add(self.textview)
        self.scrolled.set_min_content_height(300)
        vbox.pack_start(self.scrolled, True, True, 0)

        # ── Barra de progresso ─────────────────────────────────────────────
        self.progress = Gtk.ProgressBar()
        self.progress.set_show_text(True)
        self.progress.set_text("Aguardando execução...")
        vbox.pack_start(self.progress, False, False, 4)

        # ── Barra de status ────────────────────────────────────────────────
        self.status_bar = Gtk.Statusbar()
        self._ctx = self.status_bar.get_context_id("main")
        self.status_bar.push(self._ctx, "Pronto.")
        vbox.pack_start(self.status_bar, False, False, 0)

    # ---------------------------------------------------------------------- #
    # Helpers de log
    # ---------------------------------------------------------------------- #

    def log(self, mensagem: str):
        """Adiciona uma linha ao log e rola para o final."""
        end_iter = self.buffer.get_end_iter()
        self.buffer.insert(end_iter, mensagem + "\n")
        mark = self.buffer.create_mark(None, self.buffer.get_end_iter(), False)
        self.textview.scroll_to_mark(mark, 0.05, True, 0.0, 1.0)

    def set_status(self, msg: str):
        """Atualiza a barra de status."""
        self.status_bar.pop(self._ctx)
        self.status_bar.push(self._ctx, msg)

    # ---------------------------------------------------------------------- #
    # Controle de UI durante execução
    # ---------------------------------------------------------------------- #

    def _set_ui_running(self, running: bool):
        """Habilita/desabilita widgets durante a execução."""
        self._executando = running
        self.btn_run.set_sensitive(not running)
        self.btn_stop.set_sensitive(running)
        self.fc_input.set_sensitive(not running)
        self.fc_output.set_sensitive(not running)
        self.combo_module.set_sensitive(not running)
        self.combo_mode.set_sensitive(not running)

    # ---------------------------------------------------------------------- #
    # Callbacks de botões
    # ---------------------------------------------------------------------- #

    def on_run_clicked(self, widget):
        input_dir   = self.fc_input.get_filename()
        output_dir  = self.fc_output.get_filename()
        irma_module = self.combo_module.get_active_text()
        seq_mode    = self.combo_mode.get_active_id() or "auto"

        # Validação dos campos
        erros = []
        if not input_dir:
            erros.append("• Diretório de entrada não selecionado.")
        elif not os.path.isdir(input_dir):
            erros.append(f"• Diretório de entrada inválido: {input_dir}")

        if not output_dir:
            erros.append("• Diretório de saída não selecionado.")

        if not irma_module:
            erros.append("• Nenhum módulo IRMA selecionado.")

        if erros:
            self._show_error_dialog("Campos obrigatórios", "\n".join(erros))
            return

        # Cria o diretório de saída se necessário
        if not os.path.isdir(output_dir):
            try:
                os.makedirs(output_dir)
                self.log(f"📁 Diretório de saída criado: {output_dir}")
            except Exception as e:
                self._show_error_dialog("Erro", f"Não foi possível criar o diretório de saída:\n{e}")
                return

        # Localiza o script shell
        script_path = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "script_influenza_gui.sh"
        )
        if not os.path.isfile(script_path):
            self._show_error_dialog(
                "Script não encontrado",
                f"O arquivo 'script_influenza_gui.sh' não foi encontrado em:\n{script_path}\n\n"
                "Certifique-se de que ambos os arquivos estão na mesma pasta."
            )
            return

        # Salva configurações para a próxima sessão
        save_config({
            "input_dir":   input_dir,
            "output_dir":  output_dir,
            "irma_module": irma_module,
            "seq_mode":    seq_mode,
        })

        # Limpa o log e inicia
        self.buffer.set_text("")
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.log(f"{'═'*60}")
        self.log(f"  MK Flu-Pipe — Iniciado em {timestamp}")
        self.log(f"{'═'*60}")
        self.log(f"  Entrada : {input_dir}")
        self.log(f"  Saída   : {output_dir}")
        self.log(f"  Módulo  : {irma_module}")
        self.log(f"  Modo    : {seq_mode}")
        self.log(f"{'═'*60}")
        self.log("")

        self._set_ui_running(True)
        self.progress.set_fraction(0.0)
        self.progress.set_text("Executando...")
        self.set_status("Pipeline em execução...")

        # Monta o comando: o 4º argumento é o modo (omitido se "auto")
        comando = ["bash", script_path, input_dir, output_dir, irma_module]
        if seq_mode and seq_mode != "auto":
            comando.append(seq_mode)

        thread = threading.Thread(
            target=self._executar_em_thread,
            args=(comando,),
            daemon=True
        )
        thread.start()

    def on_stop_clicked(self, widget):
        if self._processo and self._processo.poll() is None:
            try:
                self._processo.terminate()
                self.log("\n⏹ Execução interrompida pelo usuário.")
            except Exception as e:
                self.log(f"⚠️  Erro ao tentar parar o processo: {e}")
        self._set_ui_running(False)
        self.progress.set_fraction(0.0)
        self.progress.set_text("Interrompido pelo usuário.")
        self.set_status("Pipeline interrompido.")

    def on_save_log_clicked(self, widget):
        dialog = Gtk.FileChooserDialog(
            title="Salvar Log como…",
            parent=self,
            action=Gtk.FileChooserAction.SAVE,
            buttons=(
                Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                Gtk.STOCK_SAVE,   Gtk.ResponseType.OK
            )
        )
        ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        dialog.set_current_name(f"log_mkflupipe_{ts}.txt")
        dialog.set_do_overwrite_confirmation(True)

        if dialog.run() == Gtk.ResponseType.OK:
            file_path = dialog.get_filename()
            start_iter = self.buffer.get_start_iter()
            end_iter   = self.buffer.get_end_iter()
            log_text   = self.buffer.get_text(start_iter, end_iter, True)
            try:
                with open(file_path, "w", encoding="utf-8") as f:
                    f.write(log_text)
                self.log(f"📄 Log salvo em: {file_path}")
                self.set_status(f"Log salvo em: {file_path}")
            except Exception as ex:
                self._show_error_dialog("Erro ao salvar", str(ex))
        dialog.destroy()

    def on_clear_log_clicked(self, widget):
        self.buffer.set_text("")
        self.set_status("Log limpo.")

    # ---------------------------------------------------------------------- #
    # Execução do script em thread separada
    # ---------------------------------------------------------------------- #

    def _executar_em_thread(self, comando: list):
        sucesso = False
        try:
            self._processo = subprocess.Popen(
                comando,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                encoding="utf-8",
                errors="replace"
            )

            for linha in self._processo.stdout:
                linha = linha.rstrip()
                GLib.idle_add(self.log, linha)
                GLib.idle_add(self.progress.pulse)

            self._processo.wait()
            sucesso = (self._processo.returncode == 0)

        except Exception as e:
            GLib.idle_add(self.log, f"❌ Erro inesperado ao executar o pipeline: {e}")

        finally:
            GLib.idle_add(self._finalizar_execucao, sucesso)

    def _finalizar_execucao(self, sucesso: bool):
        """Chamado na thread principal após o término do pipeline."""
        self._set_ui_running(False)

        if sucesso:
            self.progress.set_fraction(1.0)
            self.progress.set_text("Concluído com sucesso!")
            self.set_status("Pipeline finalizado com sucesso.")
            self._show_completion_dialog(sucesso=True)
        else:
            self.progress.set_fraction(0.0)
            self.progress.set_text("Finalizado com erros.")
            self.set_status("Pipeline finalizado com erros. Verifique o log.")
            self._show_completion_dialog(sucesso=False)

    # ---------------------------------------------------------------------- #
    # Diálogos
    # ---------------------------------------------------------------------- #

    def _show_completion_dialog(self, sucesso: bool):
        if sucesso:
            msg_type = Gtk.MessageType.INFO
            titulo   = "Pipeline concluído!"
            subtexto = "A execução foi finalizada com sucesso.\nVerifique os resultados no diretório de saída."
        else:
            msg_type = Gtk.MessageType.WARNING
            titulo   = "Pipeline finalizado com erros"
            subtexto = (
                "Ocorreram erros durante a execução.\n"
                "Verifique o log acima e o arquivo 'run_log.txt' no diretório de saída."
            )

        dialog = Gtk.MessageDialog(
            parent=self,
            flags=0,
            message_type=msg_type,
            buttons=Gtk.ButtonsType.OK,
            text=titulo
        )
        dialog.format_secondary_text(subtexto)
        dialog.run()
        dialog.destroy()

    def _show_error_dialog(self, titulo: str, mensagem: str):
        dialog = Gtk.MessageDialog(
            parent=self,
            flags=0,
            message_type=Gtk.MessageType.ERROR,
            buttons=Gtk.ButtonsType.OK,
            text=titulo
        )
        dialog.format_secondary_text(mensagem)
        dialog.run()
        dialog.destroy()

    # ---------------------------------------------------------------------- #
    # Fechamento da janela
    # ---------------------------------------------------------------------- #

    def _on_window_delete(self, widget, event):
        if self._executando and self._processo and self._processo.poll() is None:
            dialog = Gtk.MessageDialog(
                parent=self,
                flags=0,
                message_type=Gtk.MessageType.QUESTION,
                buttons=Gtk.ButtonsType.YES_NO,
                text="Pipeline em execução"
            )
            dialog.format_secondary_text(
                "O pipeline ainda está em execução.\n"
                "Deseja realmente fechar e interromper a execução?"
            )
            response = dialog.run()
            dialog.destroy()
            if response == Gtk.ResponseType.YES:
                try:
                    self._processo.terminate()
                except Exception:
                    pass
                Gtk.main_quit()
            return True  # bloqueia o fechamento se o usuário escolheu Não

        Gtk.main_quit()
        return False


# --------------------------------------------------------------------------- #
# Ponto de entrada
# --------------------------------------------------------------------------- #

def main():
    app = JanelaPrincipal()
    app.show_all()
    Gtk.main()


if __name__ == "__main__":
    main()
