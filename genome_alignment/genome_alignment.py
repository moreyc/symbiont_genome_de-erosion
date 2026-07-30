import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import sys

GENOMES = {
	'GCA_000015105.1': 'genomes/GCA_000015105.1_ASM1510v1_genomic.fasta',
	'GCA_030142985.1': 'genomes/GCA_030142985.1_ASM3014298v1_genomic.fasta',
	'GCA_021730105.1': 'genomes/GCA_021730105.1_ASM2173010v1_genomic.fasta',
	'GCA_016592615.1': 'genomes/GCA_016592615.1_ASM1659261v1_genomic.fasta',
	'GCA_016592595.1': 'genomes/GCA_016592595.1_ASM1659259v1_genomic.fasta',
	'GCA_021730205.1': 'genomes/GCA_021730205.1_ASM2173020v1_genomic.fasta',
	'GCA_021730185.1': 'genomes/GCA_021730185.1_ASM2173018v1_genomic.fasta',
	'GCA_013342985.1': 'genomes/GCA_013342985.1_ASM1334298v1_genomic.fasta',
	'GCA_021730225': 'genomes/GCA_021730225.1_ASM2173022v1_genomic.fasta',
	'GCA_021730245': 'genomes/GCA_021730245.1_ASM2173024v1_genomic.fasta',
	'GCA_000010405': 'genomes/GCA_000010405.1_ASM1040v1_genomic.fasta',
	'GCA_021730145': 'genomes/GCA_021730145.1_ASM2173014v1_genomic.fasta',
	'GCA_021730125': 'genomes/GCA_021730125.1_ASM2173012v1_genomic.fasta',
	'GCA_014706615': 'genomes/GCA_014706615.1_ASM1470661v1_genomic.fasta',
	'GCA_021730165': 'genomes/GCA_021730165.1_ASM2173016v1_genomic.fasta',
	'GCA_021730265': 'genomes/GCA_021730265.1_ASM2173026v1_genomic.fasta',
}


NAMES = {
	'GCA_000015105.1': 'Calyptogena magnifica SOB',
	'GCA_030142985.1': 'Thyasira sp. Haima SOB',
	'GCA_021730105.1': 'Pliocardia ponderosa SOB',
	'GCA_016592615.1': 'Isorropodon fossajaponicum SOB',
	'GCA_016592595.1': 'Abyssogena phaseoliformis SOB',
	'GCA_021730205.1': 'Calyptogena pacifica SOB',
	'GCA_021730185.1': 'Calyptogena rectimargo SOB',
	'GCA_013342985.1': 'Calyptogena fausta SOB',
	'GCA_021730225': 'Archivesica gigas SOB 15',
	'GCA_021730245': 'Archivesica gigas SOB 16',
	'GCA_000010405': 'Calyptogena okutanii SOB',
	'GCA_021730145': 'Archivesica soyoae SOB 4',
	'GCA_021730125': 'Archivesica soyoae SOB 5',
	'GCA_014706615': 'Archivessica marissinica SOB',
	'GCA_021730165': 'Archivesica extenta SOB',
	'GCA_021730265': 'Archivesica diagonalis SOB',
}

OUTPUT_DIR = 'mummer_output'
COORDS_DIR = f'{OUTPUT_DIR}/coords'

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(COORDS_DIR, exist_ok=True)

# MUMmer parameters
MIN_CLUSTER = 65      # Minimum length of a cluster of matches

# ==========================
# REORIENT GENOMES
# ==========================

print("\nReorienting genomes using dnaapler...")

os.makedirs('reoriented_genomes', exist_ok=True)

for i, path in GENOMES.items():
	print(f"\nReorienting {path}...")

	basename = os.path.basename(path).replace(".fasta", "")
	tmpdir = f'reoriented_genomes/{basename}.tmp'
	tmp_output = f'{tmpdir}/{basename}_reoriented.fasta'
	final_output = f'reoriented_genomes/{basename}_reoriented.fasta'

	dnaapler_cmd = [
		'dnaapler', 'all',
		'-i', path,
		'-o', tmpdir,
		'-p', basename,
		'-f'
	]

	result = subprocess.run(dnaapler_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
	if result.returncode != 0:
		print(result.stderr)
		sys.exit(1)

	os.rename(tmp_output, final_output)


# ========================== 
# RUN MUMMER ALIGNMENTS
# ==========================

species_names = list(GENOMES.keys())
species_names.reverse()  # Reverse to match row order
n = len(species_names)

print(f"\nRunning {n}x{n} = {n*n} pairwise alignments...")

for i, ref_name in enumerate(species_names):
		for j, qry_name in enumerate(species_names):
			if i < j:  # skip lower triangle
				continue
			ref_path = GENOMES[ref_name].replace(".fasta", "_reoriented.fasta").replace("genomes/", "reoriented_genomes/")
			qry_path = GENOMES[qry_name].replace(".fasta", "_reoriented.fasta").replace("genomes/", "reoriented_genomes/")
			
			prefix = f'{COORDS_DIR}/{i}_{j}'
			
			print(f"\nAligning [{i},{j}]: {ref_name} vs {qry_name}")
			
			cmd_nucmer = [
				'nucmer',
				'--prefix', prefix,
				'--maxmatch',
				ref_path,
				qry_path
			]
			
			result = subprocess.run(cmd_nucmer, 
								capture_output=True, 
								text=True)
			if result.returncode != 0:
				print(result.stderr)
				sys.exit(1)
			
			# Convert delta to coords format
			cmd_coords = [
				'show-coords',
				'-r',     # Sort by reference
				'-c',     # Show percent coverage
				'-l',     # Show alignment length
				f'{prefix}.delta'
			]
			
			with open(f'{prefix}.coords', 'w') as f:
				result = subprocess.run(cmd_coords, stdout=f, 
										stderr=subprocess.PIPE, text=True)
				if result.returncode != 0:
					print(result.stderr)
					sys.exit(1)
			print(f"\tCreated {prefix}.coords")


# ==========================
# PARSE COORDS FILES
# ==========================

def parse_coords_file(coords_file):
	"""
	Parse MUMmer show-coords output.
	
	Format:
	[S1] [E1] | [S2] [E2] | [LEN 1] [LEN 2] | [% IDY] | [LEN R] [LEN Q] | [COV R] [COV Q] | [TAGS]
	
	Returns: list of alignment blocks with coordinates and identity
	"""
	alignments = []
	
	with open(coords_file, 'r') as f:
		lines = f.readlines()

	# Skip header lines - look for the separator line with '='
		start_idx = 0
		for i, line in enumerate(lines):
			if '=' in line:
				start_idx = i + 1  # Start after the separator
				break
		
		# Parse data lines
		for line in lines[start_idx:]:
			line = line.strip()
			if not line or line.startswith('=') or line.startswith('NUCMER'):
				continue
			
			# Remove "|" characters and split
			clean_line = line.replace('|', ' ')
			parts = clean_line.split()
			
			# Expected format after cleaning:
			# [0]S1 [1]E1 [2]S2 [3]E2 [4]LEN1 [5]LEN2 [6]%IDY [7]LENR [8]LENQ [9]COVR [10]COVQ [11]TAG1 [12]TAG2

			ref_start = int(parts[0])
			ref_end = int(parts[1])
			qry_start = int(parts[2])
			qry_end = int(parts[3])
			identity = float(parts[6])
			
			alignments.append({
				'ref_start': ref_start,
				'ref_end': ref_end,
				'qry_start': qry_start,
				'qry_end': qry_end,
				'identity': identity,
				'length': abs(ref_end - ref_start)
			})
		
		return alignments

def calculate_alignment_stats(alignments):
	"""Calculate overall statistics from alignment blocks."""
	total_aligned = sum(aln['length'] for aln in alignments)
	weighted_identity = sum(aln['identity'] * aln['length'] for aln in alignments)
	
	avg_identity = weighted_identity / total_aligned if total_aligned > 0 else 0
	
	return {
		'identity': avg_identity,
		'num_blocks': len(alignments),
		'total_aligned': total_aligned
	}


# =========================
# VISUALIZE ALIGNMENTS (MEGADOTPLOT)
# =========================

n = len(species_names)

# Create figure
fig = plt.figure(figsize=(16, 16))
gs = GridSpec(n + 2, n + 1, figure=fig, 
				left=0.1, right=0.98, top=0.95, bottom=0.08,
				hspace=0.05, wspace=0.05,
				width_ratios=[1]*n + [0.05],
				height_ratios=[0.05] + [1]*n + [0.08])

# Store data for all cells
all_stats = np.zeros((n, n))

# Plot each cell
for i in range(n):
	for j in range(n):
		if i < j:  # skip lower triangle
			continue
		coords_file = f'{COORDS_DIR}/{i}_{j}.coords'
		
		# Parse alignments
		alignments = parse_coords_file(coords_file)
		stats = calculate_alignment_stats(alignments)
		all_stats[i, j] = stats['identity']
		
		# Create subplot
		ax = fig.add_subplot(gs[i + 1, j])
		
		# Plot alignment blocks
		if alignments:
			for aln in alignments:
				# Determine if reverse complement
				if aln['qry_start'] < aln['qry_end']:
					# Forward alignment
					color = plt.cm.RdYlGn(aln['identity'] / 100)
				else:
					# Reverse alignment
					color = plt.cm.Blues(aln['identity'] / 100)
				
				ax.plot([aln['ref_start'], aln['ref_end']], 
						[aln['qry_start'], aln['qry_end']], 
						color=color, linewidth=1.5, alpha=0.7)
		
		# Formatting
		ax.set_xlim(0, None)
		ax.set_ylim(0, None)
		
		if i == j:
			# Self-alignment - add gray background
			ax.set_facecolor('#f5f5f5')
		
		# Remove tick labels
		ax.set_xticklabels([])
		ax.set_xticks([])
		ax.set_yticklabels([])
		ax.set_yticks([])
		
		# Add grid
		ax.grid(True, alpha=0.2)
		
		# # Add identity text for non-self alignments
		# if i != j and stats['identity'] > 0:
		# 	ax.text(0.95, 0.05, f"{stats['identity']:.1f}%", 
		# 			transform=ax.transAxes,
		# 			ha='right', va='bottom',
		# 			fontsize=7, 
		# 			bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

# Add row labels (y-axis - query genomes)
for i, name in enumerate(NAMES.get(species_names[i], species_names[i]) for i in range(n)):
	ax = fig.add_subplot(gs[i + 1, 0])
	ax.text(-0.25, 0.5, name, rotation=0, ha='right', va='center',
			fontsize=30)
	ax.axis('off')

# Add column labels (bottom row)
for j, name in enumerate(NAMES.get(species_names[j], species_names[j]) for j in range(n)):
	ax = fig.add_subplot(gs[gs.nrows - 1, j])  # move to last row
	ax.text(0.5, 0.5, name, rotation=45,
			ha='right', va='top',
			fontsize=30)
	ax.axis('off')

# Main title
fig.suptitle('Ruthia (clade II) Genome Synteny Matrix', 
			fontsize=16, fontweight='bold', y=1)

plt.savefig(f'{OUTPUT_DIR}/alignment_matrix.svg', dpi=300, bbox_inches='tight', transparent=True)
plt.savefig(f'{OUTPUT_DIR}/alignment_matrix.pdf', bbox_inches='tight', transparent=True)
print(f"\nSaved matrix to {OUTPUT_DIR}/alignment_matrix.svg and .pdf")

