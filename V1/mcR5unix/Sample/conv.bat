# Shell script for the convolution program "conv"
# Feb. 2, 1992
#
# Format: conv.bat filename(s) output_type
# output_type includes: Rr, Ra, Az ...

# Check parm, echo the help if something is wrong.
if ($#argv == 0 || $#argv >= 3) then
  echo 'Usage: conv.bat "input_filename(s)" output_type'
  echo "output_type includes: "
  echo "I,  3,  K"
  echo "Al, Az, Arz"
  echo "Rr, Ra, Rra"
  echo "Tr, Ta, Tra"
  exit
endif

# Check the second parameter.
if (! ($2 =~ [Ii3Kk] || \
    $2 =~ [Aa][LlZz] || \
    $2 =~ [Aa][Rr][Zz] || \
    $2 =~ [RrTt][RrAa] || \
    $2 =~ [RrTt][Rr][Aa])) then
  echo "Wrong parm -- $2"
  echo "output_type includes: "
  echo "I,  3,  K"
  echo "Al, Az, Arz"
  echo "Rr, Ra, Rra"
  echo "Tr, Ta, Tra"
  exit(3)
endif

foreach infile ($1)
  # make sure the file is existent and readable.
  if (! -e $infile) then
    echo "File $infile not exist"
    exit(2);
  else if (! -r $infile) then
    echo "File $infile not readable"
    exit(2)
  endif

  # remove the existent output files, if any.
  if ( -e $infile:r.$2) then
    rm $infile:r.$2 
  endif

  # echo the command sequence to conv.
  (echo i;echo $infile;\
    echo oo;echo $2;echo $infile:r.$2;echo q;echo q;echo y)\
    |conv>/dev/null 

end # of foreach
