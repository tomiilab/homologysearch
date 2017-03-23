#!/usr/bin/perl
use strict;
use warnings;

MAIN:
{
	my $scop_q=$ARGV[0]; # ssearch_miqs
	my $scop_r=$ARGV[1]; # scop20_test
	my $standard=$ARGV[2]; # jg sf fold family
	evaluate($standard, $scop_q,$scop_r);
}


sub evaluate
{
	my $standard=$_[0];
	my $scop_q=$_[1];
	my $scop_r=$_[2];
	
	my $flag;
	if($standard eq "jg")
	{
		$flag=&Criteria($scop_q,$scop_r);
	}
	elsif($standard eq "sf" or $standard eq "fold" or $standard eq "family")
	{
		my @scop_q=split(/\./,$scop_q);
		my @scop_r=split(/\./,$scop_r);
		my $fold_q=$scop_q[0].".".$scop_q[1];
		my $fold_r=$scop_r[0].".".$scop_r[1];
		
		if($fold_q ne $fold_r)
		{
			$flag="F";
		}
		else
		{
			if($standard eq "sf")
			{
				my $sf_q=$scop_q[0].".".$scop_q[1].".".$scop_q[2];
				my $sf_r=$scop_r[0].".".$scop_r[1].".".$scop_r[2];
				if($sf_q eq $sf_r)
				{
					$flag="A";
				}
				else
				{
					$flag="B";
				}
			}
			elsif($standard eq "fold")
			{
				$flag="A";
			}
			elsif($standard eq "family")
			{
				my $family_q=$scop_q;
				my $family_r=$scop_r;
				if($family_q eq $family_r)
				{
					$flag="A";
				}
				else
				{
					$flag="B";
				}
			}
		}
	}
	print $flag;
	return $flag;

}

sub Criteria {
	#1.69 version
    #Takes as an input two SCOP classifications, and returns a flag.
    #A if they're the same, B if it's ambiguous, and F if they're different
    my $flag = "F";
    my $one  = $_[0];
    my $two  = $_[1];
    my $cf1;
    my $cf2;
    my $sf1;
    my $sf2;
    my $fa1;
    my $fa2;
    my %rossmann = (
        'c.2',  0, 'c.3',  0, 'c.4',  0, 'c.5',  0,
        'c.27', 0, 'c.28', 0, 'c.30', 0, 'c.31', 0
    );

    #these all have notes in SCOP

    if ( $one =~ /^(\w\.\d+)(\.\d+)(\.\d+)/ ) {
        $fa1 = "$1$2$3";
        $sf1 = "$1$2";
        $cf1 = "$1";
    }
    else {
        print STDERR "Error parsing classification: $one\n";
    }
    if ( $two =~ /^(\w\.\d+)(\.\d+)(\.\d+)/ ) {
        $fa2 = "$1$2$3";
        $sf2 = "$1$2";
        $cf2 = "$1";
    }
    else {
        print STDERR "Error parsing classification: $two\n";
    }

	
    #Same fold ambiguous
    if ( $cf1 eq $cf2 ) {
        $flag = "B";
    }

    #plain right
    if ( $sf1 eq $sf2 ) {
        $flag = "A";
    }

    #Fixed as of SCOP 1.63!
##Unless Membrane all-alpha
    #if (($sf1 eq 'f.2.1' or $sf2 eq 'f.2.1') and $fa1 ne $fa2){
    #$flag=F;
    #}

    #TIM barrels
    #first 7 similar -note in SCOP-
    if (
        (
               $sf1 eq 'c.1.1'
            or $sf1 eq 'c.1.2'
            or $sf1 eq 'c.1.3'
            or $sf1 eq 'c.1.4'
            or $sf1 eq 'c.1.5'
            or $sf1 eq 'c.1.6'
            or $sf1 eq 'c.1.7'
        )
        and (  $sf2 eq 'c.1.1'
            or $sf2 eq 'c.1.2'
            or $sf2 eq 'c.1.3'
            or $sf2 eq 'c.1.4'
            or $sf2 eq 'c.1.5'
            or $sf2 eq 'c.1.6'
            or $sf2 eq 'c.1.7' )
      )
    {
        $flag = "A";
    }
    elsif ( ( $cf1 eq 'c.1' and $cf2 eq 'c.1' and ( $sf1 ne $sf2 ) ) ) {
        $flag = "B";
    }

    #Rossmann-like
    if ( exists( $rossmann{$cf1} ) and exists( $rossmann{$cf2} ) ) {
        if ( $cf1 eq $cf2 ) {
            $flag = "A";
        }
        else {
            $flag = "B";
        }
    }

    #Nah! Julian thinks now it's a grower, since 1.63
##as of 1.57  c.23.12 looks like superposes OK
#if ((exists($rossmann{$cf1}) and $sf2 eq 'c.23.12') or (exists($rossmann{$cf2}) and $sf1 eq 'c.23.12')){
#$flag=1;
#}
#Old note -correspondance checked-
#1ykf, residues 151-314 superpose with 1eiz: 2.26382 ANGSTROMS/ATOM over 98 residues
    if (   ( exists( $rossmann{$cf1} ) and $cf2 eq 'c.66' )
        or ( exists( $rossmann{$cf2} ) and $cf1 eq 'c.66' ) )
    {
        $flag = "B";
    }

#1lvh superposes with 1ek6 to 2.19488 ANGSTROMS/ATOM over 88 residues(same topology), BUT is VERY different
    if (   ( exists( $rossmann{$cf1} ) and $cf2 eq 'c.108' )
        or ( exists( $rossmann{$cf2} ) and $cf1 eq 'c.108' ) )
    {
        $flag = "B";
    }

    #Can't find it and there are no cross-hits anyway.
##Old note -correspondance checked-
#if ((exists($rossmann{$cf1}) and $cf2 eq 'c.32') or (exists($rossmann{$cf2}) and $cf1 eq 'c.32')){
#$flag="B";
#}
#note: the ATP nucleotide-binding site is similar to that of the NAD-binding Rossmann-folds
    if (   ( exists( $rossmann{$cf1} ) and $sf2 eq 'c.111.1' )
        or ( exists( $rossmann{$cf2} ) and $sf1 eq 'c.111.1' ) )
    {
        $flag = "B";
    }

    #Other rules
    #beta propellors 4-8 blades
    if (
        (
               $cf1 eq 'b.66'
            or $cf1 eq 'b.67'
            or $cf1 eq 'b.68'
            or $cf1 eq 'b.69'
            or $cf1 eq 'b.70'
        )
        and (  $cf2 eq 'b.66'
            or $cf2 eq 'b.67'
            or $cf2 eq 'b.68'
            or $cf2 eq 'b.69'
            or $cf2 eq 'b.70' )
      )
    {
        if ( $cf1 eq $cf2 ) {
            $flag = "A";
        }
        else {
            $flag = "B";
        }
    }

    #Note in SCOP, Similar in architecture but partly differs in topology
    if (   ( $cf1 eq 'c.94' and $cf2 eq 'c.93' )
        or ( $cf2 eq 'c.94' and $cf1 eq 'c.93' ) )
    {
        $flag = "B";
    }

#This re-classified in 1.67 (Cutinase-like)
#if (($sf1 eq 'c.23.9' and $sf2 eq 'c.69.1') or ($sf2 eq 'c.23.9' and $sf1 eq 'c.69.1')){
#$flag="B";
#}
#This re-classified in 1.63
#if (($fa1 eq 'f.2.1.10' and $sf2 eq 'c.108.1') or ($fa2 eq 'f.2.1.10' and $sf1 eq 'c.108.1')){
#$flag="A";
#}
#Very similar alpha super-helix
    if (   ( $sf1 eq 'a.118.8' and $sf2 eq 'a.118.6' )
        or ( $sf2 eq 'a.118.8' and $sf1 eq 'a.118.6' ) )
    {
        $flag = "B";
    }

    #Similar motif sulphur binding
    if (   ( $sf1 eq 'd.58.1' and $sf2 eq 'a.1.2' )
        or ( $sf2 eq 'd.58.1' and $sf1 eq 'a.1.2' ) )
    {
        $flag = "B";
    }

    #OK note in SCOP one of the previous cases
    if (   ( $sf1 eq 'a.137.4' and $sf2 eq 'c.96.1' )
        or ( $sf2 eq 'a.137.4' and $sf1 eq 'c.96.1' ) )
    {
        $flag = "B";
    }

    #OK same fold and general look
    if (   ( $sf1 eq 'b.42.5' and $sf2 eq 'b.42.1' )
        or ( $sf2 eq 'b.42.5' and $sf1 eq 'b.42.1' ) )
    {
        $flag = "B";
    }

  #Note in SCOP. Leucine rich repeats both of them,  structures look the same OK
    if (   ( $sf1 eq 'c.10.1' and $sf2 eq 'c.10.2' )
        or ( $sf2 eq 'c.10.1' and $sf1 eq 'c.10.2' ) )
    {
        $flag = "B";
    }
##re-classified in SCOP 1.67
##Obvious sequence homology with blast,  one is beta-beta-alpha superhelix,  and one is beta-alpha togethor in PFAM
#if (($sf1 eq 'c.11.1' and $sf2 eq 'c.10.2') or ($sf2 eq 'c.11.1' and $sf1 eq 'c.10.2')){
#$flag="B";
#}
##Obvious sequence homology with blast,  one is beta-beta-alpha superhelix,  and one is beta-alpha togethor in PFAM
#if (($sf1 eq 'c.11.1' and $sf2 eq 'c.10.1') or( $sf2 eq 'c.11.1' and $sf1 eq 'c.10.1')){
#$flag="B";
#}
#Note in SCOP, contains P-loop
    if (   ( $sf1 eq 'c.91.1' and $sf2 eq 'c.37.1' )
        or ( $sf2 eq 'c.91.1' and $sf1 eq 'c.37.1' ) )
    {
        $flag = "B";
    }

    #Note in SCOP, shared motif
    if (   ( $sf1 eq 'd.51.1' and $sf2 eq 'd.52.3' )
        or ( $sf2 eq 'd.51.1' and $sf1 eq 'd.52.3' ) )
    {
        $flag = "B";
    }

    #Note in SCOP: variant of beta/alpha barrel
    if (   ( $cf1 eq 'c.6' and $cf2 eq 'c.1' )
        or ( $cf2 eq 'c.6' and $cf1 eq 'c.1' ) )
    {
        $flag = "B";
    }

    #Note in SCOP, possible link
    if (   ( $sf1 eq 'a.24.1' and $sf2 eq 'a.63.1' )
        or ( $sf2 eq 'a.24.1' and $sf1 eq 'a.63.1' ) )
    {
        $flag = "B";
    }

###NEw to 1.69
 #Shared helix no note in SCOP, sequence identical, structure different position
    if (   ( $sf1 eq 'a.8.4' and $sf2 eq 'b.130.1' )
        or ( $sf2 eq 'a.8.4' and $sf1 eq 'b.130.1' ) )
    {
        $flag = "B";
    }

#Shared helix no note in SCOP, sequence aligns, structure looks same (domain boundary definition problem)
    if (   ( $sf1 eq 'a.5.10' and $sf2 eq 'd.79.4' )
        or ( $sf2 eq 'a.5.10' and $sf1 eq 'd.79.4' ) )
    {
        $flag = "B";
    }

#Shared helix no note in SCOP, sequence alignment very bad, but includes a hyper-variable linker region. Looks like a misclassified? part of a strange chain
    if (   ( $sf1 eq 'f.21.1' and $sf2 eq 'f.32.1' )
        or ( $sf2 eq 'f.21.1' and $sf1 eq 'f.32.1' ) )
    {
        $flag = "B";
    }
###-----------

    #-----------------------------------
    return ($flag);
}
