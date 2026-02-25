import streamlit as st
import subprocess
import os
from io import StringIO, BytesIO  # ← добавлен BytesIO для Excel
import pandas as pd
import shutil
import time
from Bio import SeqIO            # ← новый импорт для конвертера

# == Константы ==
DIAMOND_DEFAULT_PATH = r"C:\Users\diamond\diamond.exe"
SUPPORTED_QUERY_TYPES = ["fasta", "fa", "faa"]
TARGET_GENES = {
    "nif_cluster": ["nifH", "nifD", "nifK", "nifE", "nifN", "nifB", "nifQ", "nifV"],
    "pho_genes":   ["phoA", "phoB", "phoD", "phoE", "phoH", "phoR", "phoU"],
    "custom": []
}
CUSTOM_TEMP_DIR = r"C:\Users\diamond_tempp"


def clear_temp_dir():
    """Очищаем временную папку и создаём заново"""
    if os.path.exists(CUSTOM_TEMP_DIR):
        try:
            shutil.rmtree(CUSTOM_TEMP_DIR)
        except Exception as e:
            st.warning(f"Не удалось очистить временную папку: {str(e)}")
    os.makedirs(CUSTOM_TEMP_DIR, exist_ok=True)
    time.sleep(2)


# ============================================================
# == Основное приложение ==
# ============================================================
def main():
    st.set_page_config(page_title="DIAMOND Suite", layout="wide", page_icon="🧬")
    
    st.sidebar.title("DIAMOND Suite")
    
    # ← НОВЫЙ режим "Genome Mining" добавлен СЮДА
    app_mode = st.sidebar.radio(
        "Select Mode",
        ["Sequence Search", "Create Database", "FNA → FAA Converter", "Genome Mining"]
    )
    
    st.sidebar.header("Configuration")
    diamond_path = st.sidebar.text_input("Path to DIAMOND executable", DIAMOND_DEFAULT_PATH)
    diamond_valid, diamond_msg = validate_diamond_path(diamond_path)
    
    if not diamond_valid:
        st.sidebar.error(diamond_msg)
        st.sidebar.warning("Try specifying full path to diamond.exe")
    else:
        st.sidebar.success(diamond_msg)
    
    # ← ОДИН блок if-elif для всех 4 режимов
    if app_mode == "Sequence Search":
        show_search_interface(diamond_path)
    elif app_mode == "Create Database":
        show_database_creator(diamond_path)
    elif app_mode == "FNA → FAA Converter":
        show_fna_to_faa_converter()
    else:  # Genome Mining
        show_genome_mining_interface(diamond_path)


# ============================================================
# == Проверка DIAMOND ==
# ============================================================
def validate_diamond_path(path):
    if not path:
        return False, "Path not specified"
    if not os.path.exists(path):
        return False, f"File not found: {path}"
    if os.path.isdir(path):
        return False, f"Path is a directory: {path}"
    try:
        result = subprocess.run([path, "version"], capture_output=True, text=True, timeout=10)
        if "diamond v" in result.stdout:
            return True, "DIAMOND is ready"
        return False, "Invalid DIAMOND executable"
    except Exception as e:
        return False, f"Check failed: {str(e)}"


# ============================================================
# == Интерфейс поиска ==
# ============================================================
def show_search_interface(diamond_path):
    st.title("🔍 DIAMOND Sequence Search")
    with st.expander("Search Settings", expanded=True):
        col1, col2 = st.columns(2)
        with col1:
            seq_type = st.radio("Sequence type", ["protein", "dna"])
            search_mode = st.radio("Gene selection mode", ["Preset genes", "Custom genes"])
            if search_mode == "Preset genes":
                gene_category = st.selectbox("Gene category", list(TARGET_GENES.keys()))
                selected_genes = st.multiselect("Select genes", TARGET_GENES[gene_category])
            else:
                custom_genes = st.text_area("Enter genes (one per line)")
                selected_genes = [g.strip() for g in custom_genes.split("\n") if g.strip()]
        with col2:
            evalue         = st.number_input("E-value threshold", min_value=0.0, value=0.001)
            max_target_seqs = st.number_input("Max target sequences", min_value=1, value=5)
            threads        = st.number_input("Threads", min_value=1, value=2, max_value=4)

    st.header("Input Data")
    col1, col2 = st.columns(2)
    with col1:
        st.subheader("Query Sequences")
        query_file = st.file_uploader("Upload query file", type=SUPPORTED_QUERY_TYPES)
        if not query_file and st.button("Generate example query"):
            example_sequence = "ATGCGTACGT" if seq_type == "dna" else "MRELEE"
            example_query = "\n".join(
                [f">gene_{gene}\n{example_sequence}" for gene in selected_genes[:3]]
            )
            st.code(example_query, language="fasta")
            query_file = StringIO(example_query)
    with col2:
        st.subheader("Database")
        db_option = st.radio("Database source", ["Local database", "Remote database (nr)"])
        db_file = st.file_uploader("Upload DIAMOND database", type=["dmnd"]) \
                  if db_option == "Local database" else None

    if st.button("🚀 Run Search", type="primary"):
        if not selected_genes:
            st.error("Please select genes to search")
        elif not query_file:
            st.error("Please upload query file")
        elif db_option == "Local database" and not db_file:
            st.error("Please upload database file")
        else:
            with st.spinner("Running search..."):
                try:
                    clear_temp_dir()
                    query_path = os.path.join(CUSTOM_TEMP_DIR, "query.fasta")
                    with open(query_path, "wb") as f:
                        if isinstance(query_file, StringIO):
                            f.write(query_file.getvalue().encode())
                        else:
                            f.write(query_file.getbuffer())

                    if db_option == "Local database":
                        db_path = os.path.join(CUSTOM_TEMP_DIR, "db.dmnd")
                        with open(db_path, "wb") as f:
                            f.write(db_file.getbuffer())
                    else:
                        db_path = "nr"

                    output_path = os.path.join(CUSTOM_TEMP_DIR, "results.tsv")
                    success, results = run_diamond_search(
                        diamond_path, query_path, db_path, output_path,
                        seq_type, evalue, max_target_seqs, threads
                    )

                    if success:
                        display_results(results)
                    else:
                        st.error(f"Search failed: {results}")

                    log_path = os.path.join(CUSTOM_TEMP_DIR, "diamond_log.txt")
                    if os.path.exists(log_path):
                        with open(log_path, "r") as f:
                            st.text_area("DIAMOND Logs", f.read(), height=200)
                except Exception as e:
                    st.error(f"Error: {str(e)}")


# ============================================================
# == Запуск DIAMOND ==
# ============================================================
def run_diamond_search(diamond_path, query_path, db_path, output_path,
                       seq_type="protein", evalue=0.001, max_target_seqs=5, threads=2):
    try:
        cmd = [
            diamond_path,
            "blastx" if seq_type == "dna" else "blastp",
            "--query", query_path,
            "--db", db_path,
            "--out", output_path,
            "--evalue", str(evalue),
            "--max-target-seqs", str(max_target_seqs),
            "--threads", str(threads),
            "--outfmt", "6",
            "--tmpdir", CUSTOM_TEMP_DIR,
            # ВАЖНО: без --log вообще
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)

        # Сохраняем stdout/stderr, чтобы можно было смотреть в интерфейсе
        log_path = os.path.join(CUSTOM_TEMP_DIR, "diamond_log.txt")
        with open(log_path, "w", encoding="utf-8") as log_file:
            log_file.write("Command:\n" + " ".join(cmd) + "\n\n")
            log_file.write(f"STDOUT:\n{result.stdout}\n\nSTDERR:\n{result.stderr}\n")

        if result.returncode != 0:
            msg = f"DIAMOND error (code {result.returncode}):\nTemporary dir: {CUSTOM_TEMP_DIR}\n"
            msg += result.stderr[:800] if result.stderr else ""
            return False, msg

        if not os.path.exists(output_path) or os.path.getsize(output_path) == 0:
            return False, f"No output file generated at {output_path}"

        columns = [
            "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore"
        ]
        results = pd.read_csv(output_path, sep="\t", header=None, names=columns)
        return True, results

    except subprocess.TimeoutExpired:
        return False, "DIAMOND timed out (300s)"
    except Exception as e:
        return False, f"Unexpected error: {str(e)}"




# ============================================================
# == Вывод результатов  ← добавлена кнопка Excel ==
# ============================================================
def display_results(results):
    st.success(f"Search completed! Found {len(results)} hits.")
    st.dataframe(results)

    col1, col2 = st.columns(2)

    with col1:
        st.download_button(
            label="📥 Download as CSV",
            data=results.to_csv(index=False).encode("utf-8"),
            file_name="results.csv",
            mime="text/csv"
        )

    with col2:
        # ← Новая кнопка: скачать Excel
        excel_buffer = BytesIO()
        with pd.ExcelWriter(excel_buffer, engine="openpyxl") as writer:
            results.to_excel(writer, index=False, sheet_name="DIAMOND Results")
        st.download_button(
            label="📊 Download as Excel",
            data=excel_buffer.getvalue(),
            file_name="results.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        )


# ============================================================
# == Создание базы DIAMOND ==
# ============================================================
def show_database_creator(diamond_path):
    st.title("🛠️ Create DIAMOND Database")
    fasta_file = st.file_uploader("Upload FASTA file", type=SUPPORTED_QUERY_TYPES)
    db_name = st.text_input("Database name", "my_database")
    if st.button("Create Database"):
        if not fasta_file:
            st.error("Please upload a FASTA file")
        else:
            with st.spinner("Creating database..."):
                try:
                    clear_temp_dir()
                    fasta_path = os.path.join(CUSTOM_TEMP_DIR, fasta_file.name)
                    with open(fasta_path, "wb") as f:
                        f.write(fasta_file.getbuffer())
                    db_path = os.path.join(CUSTOM_TEMP_DIR, f"{db_name}.dmnd")
                    result = subprocess.run(
                        [diamond_path, "makedb", "--in", fasta_path,
                         "--db", db_path, "--tmpdir", CUSTOM_TEMP_DIR],
                        capture_output=True, text=True
                    )
                    if result.returncode == 0:
                        with open(db_path, "rb") as f:
                            db_data = f.read()
                        st.success("Database created!")
                        st.download_button("Download Database", db_data, f"{db_name}.dmnd")
                    else:
                        st.error(f"Database creation failed: {result.stderr}")
                except Exception as e:
                    st.error(f"Error: {str(e)}")


# ============================================================
# == НОВАЯ ФУНКЦИЯ: FNA / FASTA → FAA Converter ==
# ============================================================
def show_fna_to_faa_converter():
    st.title("🔄 FNA / FASTA → FAA Converter")
    st.info(
        "Загрузите нуклеотидный файл (.fna / .fasta). "
        "Если файл содержит готовые CDS — выберите режим CDS. "
        "Если это геномные контиги — выберите режим Full Translation."
    )

    uploaded_file = st.file_uploader(
        "Upload nucleotide FASTA file",
        type=["fna", "fasta", "fa", "txt"]
    )

    col1, col2, col3 = st.columns(3)

    with col1:
        reading_frame = st.selectbox(
            "Reading frame",
            options=[1, 2, 3],
            index=0
        )
    with col2:
        table = st.selectbox(
            "Genetic code table",
            options=[1, 2, 4, 11],
            format_func=lambda x: {
                1:  "Standard (1)",
                2:  "Vertebrate Mitochondrial (2)",
                4:  "Mold / Protozoan / Mycoplasma (4)",
                11: "Bacterial / Archaeal / Plant Plastid (11)"
            }[x],
            index=3
        )
    with col3:
        mode = st.selectbox(
            "Translation mode",
            options=["CDS mode", "Full translation"]
        )

    min_protein_len = st.number_input(
        "Minimum protein length (aa)",
        min_value=1,
        value=10
    )

    if uploaded_file and st.button("🔄 Convert to FAA", type="primary"):
        with st.spinner("Translating sequences..."):
            try:
                content = uploaded_file.read().decode("utf-8")
                records = list(SeqIO.parse(StringIO(content), "fasta"))

                if not records:
                    st.error("Последовательности не найдены. Проверьте формат файла.")
                    return

                to_stop = (mode == "CDS mode")

                faa_lines = []
                errors = []
                skipped_short = 0

                for rec in records:
                    try:
                        seq = rec.seq[reading_frame - 1:]
                        seq = seq[: len(seq) - len(seq) % 3]

                        protein = seq.translate(
                            table=table,
                            to_stop=to_stop,
                            stop_symbol="*"
                        )

                        if len(protein) < min_protein_len:
                            skipped_short += 1
                            continue

                        desc = rec.description.replace(rec.id, "").strip()
                        header = f">{rec.id}" + (f" {desc}" if desc else "")
                        faa_lines.append(header)
                        
                        # Разбивка на строки по 60 символов
                        prot_str = str(protein)
                        faa_lines.extend(
                            prot_str[i:i+60] for i in range(0, len(prot_str), 60)
                        )
                        faa_lines.append("")  # пустая строка между записями

                    except Exception as e:
                        errors.append(f"{rec.id}: {e}")

                faa_content = "\n".join(faa_lines) + "\n"
                translated = len(records) - len(errors) - skipped_short

                # Статистика
                col_s1, col_s2, col_s3 = st.columns(3)
                col_s1.metric("Переведено", translated)
                col_s2.metric("Пропущено (короткие)", skipped_short)
                col_s3.metric("Ошибки", len(errors))

                if errors:
                    with st.expander(f"⚠️ Ошибки при трансляции ({len(errors)})"):
                        for err in errors:
                            st.warning(err)

                if translated == 0:
                    st.error(
                        "Нет результатов. Попробуйте:\n"
                        "- Уменьшить 'Minimum protein length'\n"
                        "- Переключить режим (CDS mode ↔ Full translation)\n"
                        "- Проверить рамку считывания"
                    )
                    return

                st.success(f"✅ Готово! Транслировано {translated} последовательностей.")

                # === ИСПРАВЛЕННЫЙ ПРЕДПРОСМОТР ===
                st.subheader("Preview (первые 100 строк):")
                
                preview_lines = faa_lines[:100]  # ← строго первые 100 строк
                st.code("\n".join(preview_lines), language="text")
                
                # ← Кнопка "Смотреть полную" только для полного перевода
                if mode == "Full translation" and len(faa_lines) > 100:
                    st.info("**Полная последовательность очень длинная — скачайте файл для просмотра**")
                
                # Кнопка скачивания
                base_name = os.path.splitext(uploaded_file.name)[0]
                st.download_button(
                    label="📥 Download .faa file",
                    data=faa_content.encode("utf-8"),
                    file_name=f"{base_name}.faa",
                    mime="text/plain"
                )

            except UnicodeDecodeError:
                st.error("Не удалось декодировать файл. Убедитесь, что кодировка UTF-8 или ASCII.")
            except Exception as e:
                st.error(f"Ошибка: {str(e)}")

def find_orfs(seq, min_len=30, table=11):
    """Находит реальные гены, игнорируя мусор между ними"""
    orfs = []
    for frame in range(3):
        for i in range(frame, len(seq)-2, 3):
            codon = seq[i:i+3]
            if codon in ["TAA", "TAG", "TGA"]:
                protein = seq[frame:i].translate(table=table)
                if len(protein) >= min_len:
                    orfs.append(protein)
                break  # начинаем с новой рамки после стопа
    return orfs

def show_genome_mining_interface(diamond_path):
    st.title("🧬 Genome Mining (raw FAA)")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Raw genome (FAA)")
        genome_faa = st.file_uploader(
            "Upload raw genome .faa (translated ORFs / proteins)",
            type=["faa", "fa", "fasta"]
        )

    with col2:
        st.subheader("Reference proteins")
        ref_proteins = st.file_uploader(
            "Upload reference proteins FASTA/FAA (genes you want to find)",
            type=["faa", "fa", "fasta"]
        )

    st.markdown("Настройки DIAMOND:")
    evalue = st.number_input("E-value threshold", min_value=0.0, value=1e-5, format="%.1e")
    max_target_seqs = st.number_input("Max target sequences", min_value=1, value=10)
    threads = st.number_input("Threads", min_value=1, value=2, max_value=4)

    if st.button("🚀 Run Genome Mining", type="primary"):
        if not genome_faa or not ref_proteins:
            st.error("Загрузи и сырой геном (.faa), и файл с референсными белками.")
            return

        with st.spinner("Building DIAMOND database and running search..."):
       
            genome_size = len(genome_faa.getvalue())
            ref_size = len(ref_proteins.getvalue())
        
            if genome_size < 1000:
                st.error("Сырой геном (.faa) слишком маленький — проверь формат FASTA")
                return
            if ref_size < 500:
                st.error("Reference белки слишком маленькие — нужно хотя бы несколько последовательностей")
                return
        
            st.info(f"Файлы OK: геном {genome_size//1000}kb, reference {ref_size//1000}kb")

            try:
                clear_temp_dir()

                # 1) сохраняем сырой FAA и делаем из него DIAMOND-базу
                genome_faa_path = os.path.join(CUSTOM_TEMP_DIR, "genome_raw.faa")
                with open(genome_faa_path, "wb") as f:
                    f.write(genome_faa.getbuffer())

                db_path = os.path.join(CUSTOM_TEMP_DIR, "genome_db.dmnd")
                makedb_result = subprocess.run(
                    [diamond_path, "makedb", "--in", genome_faa_path,
                     "--db", db_path, "--tmpdir", CUSTOM_TEMP_DIR],
                    capture_output=True, text=True
                )
                if makedb_result.returncode != 0:
                    st.error(f"Database creation failed:\n{makedb_result.stderr[:500]}")
                    return

                # 2) сохраняем референсные белки как query
                ref_query_path = os.path.join(CUSTOM_TEMP_DIR, "ref_proteins.faa")
                with open(ref_query_path, "wb") as f:
                    f.write(ref_proteins.getbuffer())

                # 3) запускаем diamond blastp
                output_path = os.path.join(CUSTOM_TEMP_DIR, "genome_mining.tsv")
                success, results = run_diamond_search(
                    diamond_path,
                    query_path=ref_query_path,
                    db_path=db_path,
                    output_path=output_path,
                    seq_type="protein",
                    evalue=evalue,
                    max_target_seqs=max_target_seqs,
                    threads=threads,
                )

                if success:
                    st.subheader("Results")
                    display_results(results)
                else:
                    st.error(f"Genome mining failed: {results}")

            except Exception as e:
                st.error(f"Error: {str(e)}")


# ============================================================
if __name__ == "__main__":
    main()
