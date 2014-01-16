% This is the script that tests is the Matlab can compile C code using mex

try
	mex test.c CFLAGS='\$CFLAGS -fopenmp' LDFLAGS='\$LDFLAGS -fopenmp'
	test
	fprintf('^^^^^^ check for OpenMP message\n\nMex is configured properly!\n');
	if ismac 
		!say Ready to rock!
	end
catch errr
	fprintf('Mex failed! Run mex test.c to see the error message itself\n');
end
