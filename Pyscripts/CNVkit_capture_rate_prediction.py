### MODEL FROM TOOL CNVkit

import pandas as pd
import numpy as np
import sys
	
def get_edge_bias(subarr, margin):
	"""Quantify the "edge effect" of the target tile and its neighbors.
	
	The result is proportional to the change in the target's coverage due to
	these edge effects, i.e. the expected loss of coverage near the target
	edges and, if there are close neighboring tiles, gain of coverage due
	to "spill over" reads from the neighbor tiles.
	
	(This is not the actual change in coverage. This is just a tribute.)
	"""
	tile_starts = subarr['start'].values
	tile_ends = subarr['end'].values
	tgt_sizes = tile_ends - tile_starts
	# Calculate coverage loss at (both edges of) each tile
	losses = edge_losses(tgt_sizes, margin)
	# Find tiled intervals within a margin (+/- bp) of the given probe
	# (excluding the probe itself), then calculate the relative coverage
	# "gain" due to the neighbors, if any
	gap_sizes = tile_starts[1:] - tile_ends[:-1]
	ok_gaps_mask = (gap_sizes < margin)
	ok_gaps = gap_sizes[ok_gaps_mask]
	left_gains = edge_gains(tgt_sizes[1:][ok_gaps_mask], ok_gaps, margin)
	right_gains = edge_gains(tgt_sizes[:-1][ok_gaps_mask], ok_gaps, margin)
	gains = np.zeros(len(subarr))
	gains[np.concatenate([[False], ok_gaps_mask])] += left_gains
	gains[np.concatenate([ok_gaps_mask, [False]])] += right_gains
	output_by_chrom.append(gains - losses)
	return np.concatenate(output_by_chrom)
	
	
def edge_losses(target_sizes, insert_size):
	"""Calculate coverage losses at the edges of baited regions.
	
	Letting i = insert size and t = target size, the proportional loss of
	coverage near the two edges of the baited region (combined) is:
	
	.. math :: i/2t
	
	If the "shoulders" extend outside the bait $(t < i), reduce by:
	
	.. math :: (i-t)^2 / 4it
	
	on each side, or (i-t)^2 / 2it total.
	"""
	losses = insert_size/(2 * target_sizes)
	# Drop the shoulder part that would extend past the bait
	small_mask = (target_sizes < insert_size)
	t_small = target_sizes[small_mask]
	losses[small_mask] -= ((insert_size - t_small)**2/(2 * insert_size * t_small))
	return losses
	
	
def edge_gains(target_sizes, gap_sizes, insert_size):
	import numpy as np
	"""Calculate coverage gain from neighboring baits' flanking reads.
	
	Letting i = insert size, t = target size, g = gap to neighboring bait,
	the gain of coverage due to a nearby bait, if g < i, is::
	
	.. math :: (i-g)^2 / 4it
	
	If the neighbor flank extends beyond the target (t+g < i), reduce by::
	
	.. math :: (i-t-g)^2 / 4it
	
	If a neighbor overlaps the target, treat it as adjacent (gap size 0).
	"""
	if not (gap_sizes <= insert_size).all():
		raise ValueError("Gaps greater than insert size:\n" +
		gap_sizes[gap_sizes > insert_size].head())
	gap_sizes = np.maximum(0, gap_sizes)
	gains = ((insert_size - gap_sizes)**2/(4 * insert_size * target_sizes))
	# Drop the flank part that extends past this baited region
	past_other_side_mask = (target_sizes + gap_sizes < insert_size)
	g_past = gap_sizes[past_other_side_mask]
	t_past = target_sizes[past_other_side_mask]
	gains[past_other_side_mask] -= ((insert_size - t_past - g_past)**2/(4 * insert_size * t_past))
	return gains


def gain_losses(subarray, margin):
	tile_starts = subarray['start'].values
	tile_ends = subarray['end'].values
	tgt_sizes = tile_ends - tile_starts
	losses = edge_losses(tgt_sizes, margin)
	
	gap_sizes = tile_starts[1:] - tile_ends[:-1]
	ok_gaps_mask = (gap_sizes < margin)
	ok_gaps = gap_sizes[ok_gaps_mask]

	left_gains = edge_gains(tgt_sizes[1:][ok_gaps_mask], ok_gaps, margin)
	right_gains = edge_gains(tgt_sizes[:-1][ok_gaps_mask], ok_gaps, margin)
	gains = np.zeros(len(subarray))
	gains[np.concatenate([[False], ok_gaps_mask])] += left_gains
	gains[np.concatenate([ok_gaps_mask, [False]])] += right_gains
	#return gains-losses
	tgt_sizes = tgt_sizes.tolist()
	gain_loss = gains-losses
	gain_loss = gain_loss.tolist()
	total_loss = sum([x*y for x,y in zip(tgt_sizes,gain_loss)])
	return total_loss

def target_len(df):
	tile_starts = df['start'].values
	tile_ends = df['end'].values
	tgt_sizes = tile_ends - tile_starts
	return(tgt_sizes.tolist())

def block_len(subarray):
	tile_starts = subarray['start'].values
	tile_ends = subarray['end'].values
	gap_sizes = tile_starts[1:] - tile_ends[:-1]
	return(gap_sizes.tolist())

if __name__ == '__main__':
	table = sys.argv[1]
	target_size = 60456963
	#MARGINS = range(150,410,10)
	MARGINS = range(110,410,10)
	df = pd.read_csv(table,sep='\t',header=None)
	df.columns = ['chrom','start','end']
	CHROMS = df.chrom.unique().tolist()
	for margin in MARGINS:
		output_by_chrom = []
		for ch in CHROMS:
			subarray = df[df.chrom==ch]
			#gain_loss = gain_losses(subarray, margin)
			#gain_loss = gain_loss.tolist()
			#output_by_chrom.extend(gain_loss)
			output_by_chrom.append(gain_losses(subarray, margin))
		print("%d\t%.2f" % (margin ,(target_size+sum(output_by_chrom))/target_size*100))
		#print("%d\t%.2f" % (margin ,(1+sum(output_by_chrom)/len(output_by_chrom))*100))
