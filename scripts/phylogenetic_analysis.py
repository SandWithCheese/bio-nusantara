import json
import pandas as pd
from Bio import Phylo, SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import MafftCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import re
from io import StringIO
import math # Untuk menghitung jarak

# --- Konfigurasi ---
MAFFT_EXECUTABLE = r"C:\Users\Shafi\Documents\Works\Informatika\Semester 6\Komdom\scripts\mafft\mafft.bat" # atau path lengkap ke mafft, misal "/usr/local/bin/mafft"
# script_dir_for_mafft = os.path.dirname(os.path.abspath(__file__))
# MAFFT_EXECUTABLE = os.path.join(script_dir_for_mafft, "mafft", "mafft.bat") # Contoh untuk Windows jika mafft ada di subdirektori

MAX_SPECIES_IN_SUBSET_TREE = 20

def sanitize_filename(name):
    if pd.isna(name) or name == "": return ""
    name_with_underscores = re.sub(r'\s+', '_', name)
    sanitized = re.sub(r'[^\w-]', '', name_with_underscores)
    return sanitized

def run_msa(seq_records_list, mafft_exe=MAFFT_EXECUTABLE):
    if not seq_records_list or len(seq_records_list) < 2:
        print("  MSA memerlukan minimal 2 sekuens.")
        return None
    temp_fasta_in = "temp_msa_input.fasta"
    SeqIO.write(seq_records_list, temp_fasta_in, "fasta")
    temp_fasta_out = "temp_msa_output.fasta"
    try:
        print(f"  Menjalankan MAFFT (menggunakan {mafft_exe}) untuk {len(seq_records_list)} sekuens...")
        mafft_cline = MafftCommandline(cmd=mafft_exe, input=temp_fasta_in, auto=True, quiet=True)
        stdout, stderr = mafft_cline()
        with open(temp_fasta_out, "w") as f_out: f_out.write(stdout)
        alignment = AlignIO.read(temp_fasta_out, "fasta")
        print("  MSA selesai.")
        return alignment
    except Exception as e:
        print(f"  Error saat menjalankan MAFFT: {e}")
        print("  Pastikan MAFFT terinstal dan path executable-nya benar (variabel MAFFT_EXECUTABLE).")
        return None
    finally:
        if os.path.exists(temp_fasta_in): os.remove(temp_fasta_in)
        if os.path.exists(temp_fasta_out): os.remove(temp_fasta_out)

def generate_dummy_newick_for_subset(species_names_std_subset):
    """Membuat string Newick dummy (ladderized) HANYA untuk subset spesies yang diberikan."""
    if not species_names_std_subset: return None
    if len(species_names_std_subset) == 1:
        return f"({species_names_std_subset[0]}:0.1);"
    if len(species_names_std_subset) < 2: # Seharusnya tidak terjadi jika kita memastikan subset > 1
        return f"({species_names_std_subset[0]}:0.1);"

    # Buat pohon ladderized sederhana
    # (sN:0.1, (sN-1:0.1, (... (s2:0.1, s1:0.1):0.1 ...):0.1):0.1);
    # Urutkan dulu untuk konsistensi (opsional, tapi bisa membuat struktur terlihat sama jika setnya sama)
    sorted_subset = sorted(list(set(species_names_std_subset))) # Pastikan unik dan terurut

    current_tree_segment = f"({sorted_subset[-1]}:0.1,{sorted_subset[-2]}:0.1):0.1"
    for i in range(len(sorted_subset) - 3, -1, -1):
        current_tree_segment = f"({sorted_subset[i]}:0.1,{current_tree_segment}):0.1"
    return current_tree_segment + ";"

def create_and_save_image(tree_to_draw, species_to_highlight_original, output_path, all_leaf_names_in_tree_std):
    """Fungsi umum untuk menggambar dan menyimpan pohon (bisa pohon besar atau subset dummy)."""
    num_terminals = len(tree_to_draw.get_terminals())
    fig_height = max(4, num_terminals * 0.20) # Sesuaikan penskalaan tinggi
    fig_width = 8 # Sesuaikan penskalaan lebar

    fig, ax = plt.subplots(1, 1, figsize=(fig_width, fig_height))
    
    try: # Coba root pohonnya, jika belum
        if hasattr(tree_to_draw, 'rooted') and not tree_to_draw.rooted:
            tree_to_draw.root_at_midpoint()
    except Exception as e_root:
        print(f"    Peringatan: Gagal melakukan root_at_midpoint: {e_root}")

    Phylo.draw(tree_to_draw, axes=ax, do_show=False,
               label_func=lambda x: x.name if x.is_terminal() else "",
               show_confidence=False, branch_labels=None)
    
    highlighted_in_plot = False
    species_to_highlight_std = species_to_highlight_original.replace(' ', '_')

    if species_to_highlight_std in all_leaf_names_in_tree_std:
        for txt_label in ax.texts:
            if txt_label.get_text() == species_to_highlight_std:
                txt_label.set_color('red')
                txt_label.set_fontweight('bold')
                txt_label.set_fontsize(txt_label.get_fontsize() * 1.1)
                highlighted_in_plot = True
                break
    
    title = f"Context: {species_to_highlight_original}"
    if not highlighted_in_plot and species_to_highlight_std not in all_leaf_names_in_tree_std :
        title += " (Target not in this subset view)"
    elif not highlighted_in_plot: # Target ada di tree tapi gagal dihighlight karena alasan lain
         title += " (Highlighting issue)"


    ax.set_title(title, fontsize=10)
    plt.tight_layout(rect=[0, 0.02, 1, 0.95]) # Beri sedikit ruang untuk judul
    
    try:
        plt.savefig(output_path)
        print(f"  Gambar pohon disimpan: {os.path.basename(output_path)}")
    except Exception as e:
        print(f"  Error menyimpan gambar pohon {os.path.basename(output_path)}: {e}")
    finally:
        plt.close(fig)

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    base_dir = os.path.abspath(os.path.join(script_dir, os.pardir))

    progress_file_path = os.path.join(base_dir, "data", "processed", "sequence_fetch_progress.json")
    output_image_dir = os.path.join(base_dir, "output", "phylo_trees")
    master_tree_newick_path = os.path.join(base_dir, "data", "processed", "master_phylo_tree.nwk")

    if not os.path.exists(output_image_dir):
        os.makedirs(output_image_dir); print(f"Direktori output dibuat: {output_image_dir}")

    if not os.path.exists(progress_file_path):
        print(f"Error: File progres {progress_file_path} tidak ditemukan."); return
    
    try:
        with open(progress_file_path, 'r') as f:
            content = f.read(); species_sequence_data = json.loads(content) if content else {}
        print(f"Berhasil memuat {len(species_sequence_data)} entri dari {os.path.basename(progress_file_path)}")
    except Exception as e:
        print(f"Error memuat file progres {progress_file_path}: {e}"); return

    if not species_sequence_data: print("Tidak ada data spesies dalam file progres."); return

    valid_sequence_stati = ["Not Found or Error", "Sequence Fetching Skipped", None, ""]
    seq_records_for_master_tree = []
    species_name_map_std_to_original = {}

    print("Mempersiapkan sekuens untuk pohon master...")
    for species_name_original, seq_data_fasta in species_sequence_data.items():
        if seq_data_fasta not in valid_sequence_stati and isinstance(seq_data_fasta, str) and len(seq_data_fasta) > 50:
            try:
                seq_string = (seq_data_fasta.split('\n', 1)[1].replace('\n', '') if seq_data_fasta.startswith(">") and len(seq_data_fasta.split('\n', 1)) > 1 else seq_data_fasta.replace('\n', ''))
                if not seq_string: continue
                cleaned_seq_string = re.sub(r'[^ATGCatgcNn-]', '', seq_string).upper()
                if len(cleaned_seq_string) < 50: continue
                species_name_std = species_name_original.replace(' ', '_')
                seq_records_for_master_tree.append(SeqRecord(Seq(cleaned_seq_string), id=species_name_std, description=""))
                species_name_map_std_to_original[species_name_std] = species_name_original
            except Exception as e_rec:
                print(f"  Error memproses sekuens untuk '{species_name_original}': {e_rec}")
    
    if len(seq_records_for_master_tree) < 2:
        print(f"Ditemukan kurang dari 2 spesies ({len(seq_records_for_master_tree)}) dengan sekuens valid untuk membuat pohon master."); return
    
    print(f"Total {len(seq_records_for_master_tree)} sekuens disiapkan untuk pohon master.")

    master_tree = None
    if os.path.exists(master_tree_newick_path):
        try:
            print(f"Memuat pohon master dari file Newick: {master_tree_newick_path}...")
            master_tree = Phylo.read(master_tree_newick_path, "newick")
            print("Pohon master berhasil dimuat.")
            # Verifikasi apakah semua spesies yang ada di seq_records_for_master_tree ada di master_tree
            # Ini penting jika file progres diperbarui setelah master_tree.nwk dibuat.
            # Untuk kesederhanaan, kita asumsikan master_tree.nwk selalu sinkron atau dibuat ulang jika perlu.
            # Jika tidak sinkron, sebaiknya master_tree.nwk dihapus agar dibuat ulang.
        except Exception as e_load_nwk:
            print(f"Error memuat pohon master dari Newick: {e_load_nwk}. Akan coba membuat ulang.")
            master_tree = None # Paksa pembuatan ulang

    if not master_tree:
        print("Membuat pohon master baru...")
        alignment = run_msa(seq_records_for_master_tree)
        if not alignment: print("MSA gagal. Tidak bisa membuat pohon master."); return
        
        print("Menghitung matriks jarak untuk pohon master...")
        try:
            calculator = DistanceCalculator('identity'); dist_matrix = calculator.get_distance(alignment)
            print("Matriks jarak pohon master berhasil dihitung.")
        except Exception as e_dist: print(f"Error menghitung matriks jarak pohon master: {e_dist}"); return

        print("Membuat pohon master menggunakan Neighbor Joining...")
        try:
            constructor = DistanceTreeConstructor(calculator, 'nj'); master_tree = constructor.build_tree(alignment)
            print("Pohon master berhasil dibuat.")
            # Simpan pohon master yang baru dibuat ke file Newick
            try:
                master_tree_dir = os.path.dirname(master_tree_newick_path)
                if master_tree_dir and not os.path.exists(master_tree_dir): os.makedirs(master_tree_dir)
                Phylo.write(master_tree, master_tree_newick_path, "newick")
                print(f"Pohon master disimpan ke: {master_tree_newick_path}")
            except Exception as e_save_nwk:
                print(f"Error menyimpan pohon master ke Newick: {e_save_nwk}")
        except Exception as e_tree: print(f"Error membuat pohon master: {e_tree}"); return

    if not master_tree: print("Gagal membuat atau memuat pohon master."); return

    all_leaves_master_tree_std = {term.name for term in master_tree.get_terminals() if term.name}

    print(f"\nMemulai pembuatan gambar pohon subset (maks {MAX_SPECIES_IN_SUBSET_TREE} spesies)...")
    # Iterasi menggunakan nama asli (kunci dari species_name_map_std_to_original, atau dari seq_records)
    list_of_original_target_species = sorted(list(species_name_map_std_to_original.values()))

    for target_species_original in list_of_original_target_species:
        target_species_std = target_species_original.replace(' ', '_')

        if target_species_std not in all_leaves_master_tree_std:
            print(f"  Spesies target '{target_species_original}' tidak ditemukan di pohon master, dilewati.")
            continue

        print(f"Memproses gambar untuk: {target_species_original}")
        
        # Hitung jarak dari target ke semua daun lain di master_tree
        distances = []
        try:
            target_node = next(n for n in master_tree.get_terminals() if n.name == target_species_std)
            for leaf in master_tree.get_terminals():
                if leaf.name and leaf.name != target_species_std :
                    # master_tree.distance terkadang butuh rooting eksplisit atau bisa error
                    # Jika error, kita bisa mengabaikan jarak dan hanya ambil N pertama secara acak/alfabetis
                    # Atau, jika pohon besar dan berakar, ini seharusnya bekerja.
                    try:
                        dist = master_tree.distance(target_node, leaf)
                        distances.append((leaf.name, dist))
                    except Exception as e_dist_leaf:
                        # print(f"    Tidak bisa menghitung jarak ke {leaf.name}: {e_dist_leaf}")
                        # Fallback jika jarak tidak bisa dihitung: beri jarak besar agar tidak terpilih
                        distances.append((leaf.name, float('inf'))) 
            
            distances.sort(key=lambda x: x[1]) # Urutkan berdasarkan jarak
            
            # Ambil MAX_SPECIES_IN_SUBSET_TREE - 1 tetangga terdekat
            closest_neighbors_std = [d[0] for d in distances[:MAX_SPECIES_IN_SUBSET_TREE - 1]]
            subset_species_std = sorted(list(set([target_species_std] + closest_neighbors_std)))

        except StopIteration: # Jika target_node tidak ditemukan (seharusnya tidak terjadi jika sudah dicek)
            print(f"    Spesies target '{target_species_std}' tidak ditemukan sebagai terminal di master_tree (aneh).")
            subset_species_std = [target_species_std] # Hanya spesies itu sendiri
        except Exception as e_calc_dist: # Error umum saat kalkulasi jarak
            print(f"    Error saat menghitung jarak untuk {target_species_std}: {e_calc_dist}. Menggunakan subset acak/terbatas.")
            # Fallback: ambil N pertama dari master tree jika error, atau hanya target
            temp_leaves = list(all_leaves_master_tree_std)
            if target_species_std in temp_leaves: temp_leaves.remove(target_species_std)
            subset_species_std = sorted(list(set([target_species_std] + temp_leaves[:MAX_SPECIES_IN_SUBSET_TREE -1 ])))


        if len(subset_species_std) < 1: # Harus ada minimal 1
             subset_species_std = [target_species_std] if target_species_std else []


        # Buat pohon dummy HANYA untuk subset ini
        subset_newick_string = generate_dummy_newick_for_subset(subset_species_std)
        subset_tree_to_draw = None
        if subset_newick_string:
            try:
                subset_tree_to_draw = Phylo.read(StringIO(subset_newick_string), "newick")
            except Exception as e_parse_subset:
                print(f"    Error parsing Newick untuk subset {target_species_original}: {e_parse_subset}")
        
        if not subset_tree_to_draw: # Jika gagal buat pohon subset, buat pohon 1 daun
            single_leaf_newick = f"({target_species_std}:0.1);"
            subset_tree_to_draw = Phylo.read(StringIO(single_leaf_newick), "newick")
            subset_species_std = [target_species_std] # Update daftar daun untuk visualisasi

        # Siapkan untuk penyimpanan gambar
        sanitized_basename = sanitize_filename(target_species_original)
        if not sanitized_basename: continue
        image_filename = f"{sanitized_basename}_phylo.png"
        output_image_path = os.path.join(output_image_dir, image_filename)
        
        create_and_save_image(subset_tree_to_draw, target_species_original, output_image_path, subset_species_std)

    print(f"\nSelesai membuat semua gambar pohon di {output_image_dir}")

if __name__ == "__main__":
    if os.environ.get('DISPLAY','') == '':
        print('Tidak ada display terdeteksi, menggunakan backend Agg untuk Matplotlib.')
        matplotlib.use('Agg')
    main()
