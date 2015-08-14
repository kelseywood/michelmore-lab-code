#!/usr/bin/perl
use strict; use warnings;
use Getopt::Std;
use vars qw($opt_d $opt_x $opt_y);
getopts('dx:y:');

my $DX = 20;
my $DY = 10;

die "
usage: procgff.pl [options] <gff>
options:
  -d       debug mode: no acutal clustering
  -x <int> external comparision distance [$DX]
  -y <int> internal comparision distance [$DY]
" unless @ARGV == 1;

my $DEBUG = $opt_d;
$DX = $opt_x if $opt_x;
$DY = $opt_y if $opt_y;


############
# Read GFF #
############
my %GFF;
my @Output;
while (<>) {
	if (/^##FASTA/) {
		last;
	} elsif (/^#/) {
		print;
		next;
	}
	
	my ($dna, $src, $fea, $beg, $end, $scr, $str, $frm, $grp, $oth)
		= split(/\t/, $_);
	if (not defined $grp) {
		print $_;
		die "not enough fields\n";
	}
	if (defined $oth) {
		print $_;
		die "too many fields\n";
	}
	
	my $thing = {
		fea => $fea,
		beg => $beg,
		end => $end,
		str => $str,
		gff => $_,
	};
	
	my @attr = split(/;/, $grp);
	my ($parent, $id);
	foreach my $attr (@attr) {
		if ($attr =~ /^ID=(\S+)/) {$id = $1}
		if ($attr =~ /^Parent=(\S+)/) {$parent = $1}
	}

	if ($src =~ /blast|2genome/ and $fea =~ /match/) {
		if (defined $parent) {
			push @{$GFF{$src}{$parent}{children}}, $thing;
		} elsif (defined $id) {
			if (defined $GFF{$src}{$id}{parent}) {
				die "already defined?";
			}
			$GFF{$src}{$id}{parent} = $thing;
		} else {
			die "wtf";
		}
	} elsif ($src eq 'RepeatProteinMask') {
	#	next; # skip these
	} elsif ($src eq 'UCD') {
		$_ =~ s/region/match_part/;
		push @Output, $_;
	}else {
		push @Output, $_;
	}
}


###################
# Cluster matches #
###################
foreach my $src (keys %GFF) {
	my @cluster = cluster($GFF{$src});
	foreach my $cluster (@cluster) {
		push @Output, $cluster->{parent}{gff};
		foreach my $child (@{$cluster->{children}}) {
			push @Output, $child->{gff};
		}
	}
}

print @Output;
print "##FASTA\n";
while (<>) {print}


##############################################################################
# subroutines
##############################################################################

sub cluster {
	my ($fs) = @_;
	
	# reorganizing a little
	my @gene;
	foreach my $id (keys %$fs) {
	
		push @gene, {
			beg => $fs->{$id}{parent}{beg},
			end => $fs->{$id}{parent}{end},
			str => $fs->{$id}{parent}{str},
			exons => scalar @{$fs->{$id}{children}},
			parent => $fs->{$id}{parent},
			children => $fs->{$id}{children},
		}
	}
	
	return @gene if $DEBUG;
		
	# clustering
	my @cluster;
	while (@gene) {
		my $node = shift @gene;
		my $cluster;
		for (my $i = 0; $i < @gene; $i++) {
			my $cmp = compare($node, $gene[$i]);
			if ($cmp eq 'super' or $cmp eq 'equal') {
				$cluster = $node; # combine($node, $gene[$i]);
				splice(@gene, $i);
				last;
			} elsif ($cmp eq 'subset') {
				$cluster = $gene[$i]; # combine($gene[$i], $node);
				splice(@gene, $i);
				last;
			}
		}
		if ($cluster) {
			unshift @gene, $cluster;
		} else {
			push @cluster, $node;
		}
	}
	
	return @cluster;
}

sub compare {
	my ($g1, $g2) = @_;
		
	# outer checks
	return 'none' if $g1->{exons} != $g2->{exons};
	return 'none' if $g1->{str}   ne $g2->{str};
	return 'none' if dissimilar($g1, $g2, $DX);
	
	# inner checks
	if ($g1->{exons} > 1) {
		my @c1 = @{$g1->{children}};
		my @c2 = @{$g2->{children}};
	
		# first and last match check only one side of feature
		my $dbeg = abs $c1[0]{end} - $c2[0]{end};
		my $dend = abs $c1[@c1-1]{beg} - $c2[@c2-1]{beg};
		return 'none' if $dbeg > $DY;
		return 'none' if $dend > $DY;
		
		# middle matches check both sides of feature
		for (my $i = 1; $i < @c1 -1; $i++) {
			return 'none' if dissimilar($c1[$i], $c2[$i], $DY);
		}
	}
	
	if ($g1->{beg} <= $g2->{beg} and $g1->{end} >= $g2->{end}) {
		return 'super';
	} elsif ($g1->{beg} > $g2->{beg} and $g1->{end} < $g2->{end}) {
		return 'subset';
	} else {
		return 'equal';
	}
}

sub dissimilar {
	my ($f1, $f2, $d) = @_;
	my $x = abs($f1->{beg} - $f2->{beg});
	my $y = abs($f2->{end} - $f2->{end});
	if ($x > $d or $y > $d) {return 1}
	else                    {return 0}
}

sub combine {
	my ($parent, $child) = @_;
	push @{$parent->{cluster}}, $child;
	return $parent;
}
