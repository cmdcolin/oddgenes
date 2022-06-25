use Data::Dumper;
my %hash=();
while (my $text=<>) {
    chomp($text);
    if(substr($text, 0, 1) ne '>') {
        foreach $char (split //, $text) {
            if($hash{$char}) {
                $hash{$char}++;
            } else {
                $hash{$char}=1
            }
        }
    }
}
print Dumper(\%hash)
