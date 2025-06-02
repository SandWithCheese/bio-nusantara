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
import math

# --- Konfigurasi ---
MAFFT_EXECUTABLE = r"C:\Users\Shafi\Documents\Works\Informatika\Semester 6\Komdom\scripts\mafft\mafft.bat" 
# Contoh untuk Windows jika mafft ada di subdirektori 'scripts/mafft':
# script_dir_for_mafft = os.path.dirname(os.path.abspath(__file__))
# MAFFT_EXECUTABLE = os.path.join(script_dir_for_mafft, "mafft", "mafft.bat")

MAX_SPECIES_IN_SUBSET_TREE = 20

def standardize_species_name_for_id(original_name):
    """
    Membersihkan dan menstandarisasi nama spesies menjadi format Genus_spesies.
    Menghilangkan penulis, tahun, dan sub-spesies untuk ID internal pohon.
    """
    if pd.isna(original_name) or not str(original_name).strip():
        return None

    name = str(original_name)
    # 1. Hapus konten dalam tanda kurung (seringkali penulis, tahun)
    name = re.sub(r'\s*\(.*\)\s*', '', name).strip()
    # 2. Hapus tahun di akhir jika formatnya ", Tahun" atau " Tahun" atau hanya "YYYY"
    name = re.sub(r'[, ]*\b\d{4}\b$', '', name).strip()
    # 3. Hapus "sp.", "spp.", "cf.", "aff." dan variasinya (termasuk titik)
    name = re.sub(r'\s*\b(spp?|cf|aff)\.?\b.*', '', name, flags=re.IGNORECASE).strip() # Hapus juga sisa string setelahnya
    # 4. Hapus karakter non-alfanumerik di akhir yang mungkin tersisa
    name = re.sub(r'[^a-zA-Z0-9\s]+$', '', name).strip()
    
    name_parts = name.split()
    
    if len(name_parts) >= 2:
        # Ambil Genus dan spesies, gabungkan dengan underscore
        # Pastikan Genus diawali huruf besar, spesies huruf kecil (umumnya)
        genus = name_parts[0].capitalize()
        species = name_parts[1].lower()
        return f"{genus}_{species}"
    elif len(name_parts) == 1:
        # Jika hanya satu kata (misalnya, Genus saja atau nama tunggal)
        return name_parts[0].capitalize() 
    else:
        # Jika nama menjadi kosong setelah pembersihan
        return None

def sanitize_filename(original_name_for_file):
    """
    Membersihkan nama spesies asli (dengan penulis/tahun) untuk nama file PNG.
    """
    if pd.isna(original_name_for_file) or not str(original_name_for_file).strip(): 
        return "unknown_species" # Default jika nama asli kosong
    
    name = str(original_name_for_file)
    # Ganti spasi dan karakter non-alfanumerik umum dengan underscore
    name = re.sub(r'\s+', '_', name)
    # Hapus tanda kurung (seringkali berisi koma yang sudah jadi underscore)
    name = re.sub(r'_\(_.*?_\)_', '_', name) # Pola seperti _(Penulis,_Tahun)_
    name = re.sub(r'\(.*?\)', '', name)     # Pola umum (apa saja)
    # Hapus karakter yang tidak valid untuk nama file, kecuali underscore dan hyphen
    sanitized = re.sub(r'[^\w.-]', '', name) # Izinkan titik juga
    sanitized = re.sub(r'_+', '_', sanitized) # Ganti multiple underscore dengan satu
    return sanitized.strip('_-.') or "sanitized_species_name"


def run_msa(seq_records_list, mafft_exe=MAFFT_EXECUTABLE):
    if not seq_records_list or len(seq_records_list) < 2:
        print("  MSA memerlukan minimal 2 sekuens.")
        return None
    temp_fasta_in = "temp_msa_input.fasta"
    
    unique_records_for_msa = []
    seen_ids_msa = set()
    for rec in seq_records_list:
        if rec.id not in seen_ids_msa:
            unique_records_for_msa.append(rec)
            seen_ids_msa.add(rec.id)
        else:
            print(f"    Peringatan MSA: ID duplikat '{rec.id}' ditemukan, hanya satu yang digunakan untuk MSA.")

    if len(unique_records_for_msa) < 2:
        print("  MSA memerlukan minimal 2 sekuens unik setelah deduplikasi ID untuk MSA.")
        return None

    SeqIO.write(unique_records_for_msa, temp_fasta_in, "fasta")
    temp_fasta_out = "temp_msa_output.fasta"
    try:
        print(f"  Menjalankan MAFFT (menggunakan {mafft_exe}) untuk {len(unique_records_for_msa)} sekuens unik...")
        mafft_cline = MafftCommandline(cmd=mafft_exe, input=temp_fasta_in, auto=True, quiet=True)
        stdout, stderr = mafft_cline()
        with open(temp_fasta_out, "w") as f_out: f_out.write(stdout)
        alignment = AlignIO.read(temp_fasta_out, "fasta")
        print("  MSA selesai.")
        return alignment
    except FileNotFoundError:
        print(f"  Error: MAFFT executable '{mafft_exe}' tidak ditemukan.")
        print("  Pastikan MAFFT terinstal dan MAFFT_EXECUTABLE diatur dengan benar di skrip.")
        return None
    except Exception as e:
        print(f"  Error saat menjalankan MAFFT: {e}")
        return None
    finally:
        if os.path.exists(temp_fasta_in): os.remove(temp_fasta_in)
        if os.path.exists(temp_fasta_out): os.remove(temp_fasta_out)

def generate_dummy_newick_for_subset(species_names_std_subset):
    if not species_names_std_subset: return None
    unique_ordered_subset = []
    seen = set()
    for name in species_names_std_subset: # Input sudah distandarisasi dan diurutkan
        if name not in seen: unique_ordered_subset.append(name); seen.add(name)
    if not unique_ordered_subset: return None
    if len(unique_ordered_subset) == 1: return f"({unique_ordered_subset[0]}:0.1);"
    # Perlu minimal 2 untuk membuat pasangan. Jika kurang dari 2 setelah unik, kembali ke 1.
    if len(unique_ordered_subset) < 2 : 
        return f"({unique_ordered_subset[0]}:0.1);" if unique_ordered_subset else None


    current_tree_segment = f"({unique_ordered_subset[-1]}:0.1,{unique_ordered_subset[-2]}:0.1):0.1"
    for i in range(len(unique_ordered_subset) - 3, -1, -1):
        current_tree_segment = f"({unique_ordered_subset[i]}:0.1,{current_tree_segment}):0.1"
    return current_tree_segment + ";"

def create_and_save_image(tree_to_draw, species_to_highlight_original, output_path, all_leaf_names_in_tree_std):
    num_terminals = len(tree_to_draw.get_terminals())
    fig_height = max(4, num_terminals * 0.25) 
    fig_width = 10

    fig, ax = plt.subplots(1, 1, figsize=(fig_width, fig_height))
    
    try:
        if hasattr(tree_to_draw, 'rooted') and not tree_to_draw.rooted:
            tree_to_draw.root_at_midpoint()
    except Exception: pass

    target_std_for_match = standardize_species_name_for_id(species_to_highlight_original)

    def label_prop_func(clade_name_from_tree): # clade_name_from_tree adalah ID standar
        if clade_name_from_tree == target_std_for_match:
            return {'color': 'red', 'weight': 'bold', 'size': plt.rcParams['font.size'] * 1.15}
        return {'color': 'black', 'weight': 'normal', 'size': plt.rcParams['font.size']}

    Phylo.draw(tree_to_draw, axes=ax, do_show=False,
               label_func=lambda x: x.name if x.is_terminal() else "",
               label_colors=lambda name_in_tree: label_prop_func(name_in_tree).get('color'),
               show_confidence=False, branch_labels=None)

    highlighted_in_plot = False
    if target_std_for_match in all_leaf_names_in_tree_std:
        for txt_label in ax.texts:
            # txt_label.get_text() adalah ID standar dari pohon
            if txt_label.get_text() == target_std_for_match:
                props = label_prop_func(txt_label.get_text())
                txt_label.set_fontweight(props.get('weight'))
                txt_label.set_fontsize(props.get('size'))
                highlighted_in_plot = True
                break 
    
    title_text = species_to_highlight_original # Tampilkan nama asli lengkap di judul
    title = f"Konteks Filogenetik untuk:\n{title_text}" # \n untuk potensi nama panjang
    if not highlighted_in_plot and target_std_for_match:
        title += "\n(Target tidak tersorot. Periksa kecocokan nama.)"
        print(f"    Peringatan: Gagal menyorot '{species_to_highlight_original}'.")
        print(f"      Nama standar dicari: '{target_std_for_match}'")
        print(f"      Daun tersedia di pohon gambar (subset): {all_leaf_names_in_tree_std[:5]}...")


    ax.set_title(title, fontsize=9, loc='left') # Perkecil font judul
    plt.subplots_adjust(left=0.1, right=0.95, top=0.88, bottom=0.05) # Sesuaikan margin, beri ruang lebih untuk judul
    
    try:
        plt.savefig(output_path, dpi=150)
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
    for species_name_original_from_json, seq_data_fasta in species_sequence_data.items():
        if seq_data_fasta not in valid_sequence_stati and isinstance(seq_data_fasta, str) and len(seq_data_fasta) > 50:
            try:
                seq_string = (seq_data_fasta.split('\n', 1)[1].replace('\n', '') if seq_data_fasta.startswith(">") and len(seq_data_fasta.split('\n', 1)) > 1 else seq_data_fasta.replace('\n', ''))
                if not seq_string: continue
                cleaned_seq_string = re.sub(r'[^ATGCatgcNn-]', '', seq_string).upper()
                if len(cleaned_seq_string) < 50: continue
                
                species_name_std = standardize_species_name_for_id(species_name_original_from_json)
                if not species_name_std:
                    # print(f"  Gagal menstandarisasi nama untuk ID: '{species_name_original_from_json}', dilewati.")
                    continue
                
                if species_name_std not in species_name_map_std_to_original:
                    seq_records_for_master_tree.append(SeqRecord(Seq(cleaned_seq_string), id=species_name_std, description=""))
                    species_name_map_std_to_original[species_name_std] = species_name_original_from_json
            except Exception as e_rec:
                print(f"  Error memproses sekuens untuk '{species_name_original_from_json}': {e_rec}")
    
    if len(seq_records_for_master_tree) < 2:
        print(f"Ditemukan kurang dari 2 spesies ({len(seq_records_for_master_tree)}) dengan sekuens valid untuk membuat pohon master."); return
    
    print(f"Total {len(seq_records_for_master_tree)} sekuens unik (berdasarkan ID standar) disiapkan untuk pohon master.")

    master_tree = None
    if os.path.exists(master_tree_newick_path):
        try:
            print(f"Memuat pohon master dari file Newick: {master_tree_newick_path}...")
            master_tree = Phylo.read(master_tree_newick_path, "newick")
            master_tree_leaves = {term.name for term in master_tree.get_terminals()}
            current_record_ids = {rec.id for rec in seq_records_for_master_tree}
            # Jika set daun berbeda, paksa regenerasi
            if master_tree_leaves != current_record_ids:
                 print("    Peringatan: Daun di master_tree.nwk berbeda dengan sekuens saat ini. Pohon master akan dibuat ulang.")
                 master_tree = None # Paksa regenerasi
            else:
                 print("Pohon master berhasil dimuat dan sinkron.")
        except Exception as e_load_nwk:
            print(f"Error memuat pohon master dari Newick: {e_load_nwk}. Akan coba membuat ulang.")
            master_tree = None

    if not master_tree:
        print("Membuat pohon master baru...")
        alignment = run_msa(seq_records_for_master_tree)
        if not alignment: print("MSA gagal. Tidak bisa membuat pohon master."); return
        
        print("Menghitung matriks jarak untuk pohon master...")
        try:
            calculator = DistanceCalculator('identity'); dist_matrix = calculator.get_distance(alignment)
        except Exception as e_dist: print(f"Error menghitung matriks jarak pohon master: {e_dist}"); return

        print("Membuat pohon master menggunakan Neighbor Joining...")
        try:
            constructor = DistanceTreeConstructor(calculator, 'nj'); master_tree = constructor.build_tree(alignment)
            print("Pohon master berhasil dibuat.")
            try:
                master_tree_dir = os.path.dirname(master_tree_newick_path)
                if master_tree_dir and not os.path.exists(master_tree_dir): os.makedirs(master_tree_dir)
                Phylo.write(master_tree, master_tree_newick_path, "newick")
                print(f"Pohon master disimpan ke: {master_tree_newick_path}")
            except Exception as e_save_nwk: print(f"Error menyimpan pohon master ke Newick: {e_save_nwk}")
        except Exception as e_tree: print(f"Error membuat pohon master: {e_tree}"); return

    if not master_tree: print("Gagal membuat atau memuat pohon master."); return

    all_leaves_master_tree_std = {term.name for term in master_tree.get_terminals() if term.name}
    list_of_target_species_std_keys = sorted(list(species_name_map_std_to_original.keys()))

    print(f"\nMemulai pembuatan gambar pohon subset (maks {MAX_SPECIES_IN_SUBSET_TREE} spesies)...")
    for target_species_std in list_of_target_species_std_keys:
        target_species_original = species_name_map_std_to_original[target_species_std]

        if target_species_std not in all_leaves_master_tree_std:
            print(f"  Spesies target '{target_species_original}' (standar: {target_species_std}) tidak ditemukan di pohon master setelah dimuat/dibuat, dilewati.")
            continue

        print(f"Memproses gambar untuk: {target_species_original} (ID standar: {target_species_std})")
        
        distances = []
        try:
            target_node = next(n for n in master_tree.get_terminals() if n.name == target_species_std)
            for leaf in master_tree.get_terminals():
                if leaf.name and leaf.name != target_species_std :
                    try:
                        if hasattr(master_tree, 'rooted') and not master_tree.rooted: master_tree.root_at_midpoint()
                        dist = master_tree.distance(target_node, leaf)
                        distances.append((leaf.name, dist))
                    except Exception: distances.append((leaf.name, float('inf'))) 
            
            distances.sort(key=lambda x: x[1])
            closest_neighbors_std = [d[0] for d in distances[:MAX_SPECIES_IN_SUBSET_TREE - 1]]
            subset_species_std_for_dummy = sorted(list(set([target_species_std] + closest_neighbors_std)))

        except Exception as e_calc_dist:
            print(f"    Error saat menghitung jarak/subset untuk {target_species_std}: {e_calc_dist}. Menggunakan subset terbatas.")
            temp_leaves = list(all_leaves_master_tree_std)
            if target_species_std in temp_leaves: temp_leaves.remove(target_species_std)
            import random; random.shuffle(temp_leaves) 
            subset_species_std_for_dummy = sorted(list(set([target_species_std] + temp_leaves[:MAX_SPECIES_IN_SUBSET_TREE -1 ])))

        if not subset_species_std_for_dummy : subset_species_std_for_dummy = [target_species_std] if target_species_std else []
        if not subset_species_std_for_dummy : 
            print(f"    Tidak ada spesies untuk subset pohon dummy target: {target_species_original}"); continue
        
        subset_newick_string = generate_dummy_newick_for_subset(subset_species_std_for_dummy)
        subset_tree_to_draw = None
        if subset_newick_string:
            try: subset_tree_to_draw = Phylo.read(StringIO(subset_newick_string), "newick")
            except Exception as e_parse_subset: print(f"    Error parsing Newick untuk subset {target_species_original}: {e_parse_subset}")
        
        if not subset_tree_to_draw:
            single_leaf_newick = f"({target_species_std}:0.1);" if target_species_std else "(unknown:0.1);"
            try: subset_tree_to_draw = Phylo.read(StringIO(single_leaf_newick), "newick")
            except: print(f"    Gagal membuat pohon bahkan untuk satu daun: {target_species_std}"); continue
            subset_species_std_for_dummy = [target_species_std] if target_species_std else ["unknown"]

        sanitized_filename_for_output = sanitize_filename(target_species_original) 
        if not sanitized_filename_for_output: print(f"  Nama file tidak valid untuk {target_species_original}"); continue
        image_filename = f"{sanitized_filename_for_output}_phylo.png"
        output_image_path = os.path.join(output_image_dir, image_filename)
        
        create_and_save_image(subset_tree_to_draw, target_species_original, output_image_path, subset_species_std_for_dummy)

    print(f"\nSelesai membuat semua gambar pohon di {output_image_dir}")

if __name__ == "__main__":
    if os.environ.get('DISPLAY','') == '':
        print('Tidak ada display terdeteksi, menggunakan backend Agg untuk Matplotlib.')
        matplotlib.use('Agg')
    main()

