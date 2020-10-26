restart

pos = toList((0,0)..(2,2));
patterns = sort subsets pos;
patterns = apply(patterns,sort);
# patterns
perm3 = permutations{0,1,2};
perm2 = permutations{0,1};

--function for permuting a zero pattern
--the inputs are the zero pattern, row permutation, column permutation and transposition
permutePattern = (pattern, rowPerm, colPerm, trans) -> (
    newPattern1 = apply(pattern,(x,y)->(rowPerm#x,colPerm#y));
    apply(newPattern1,a -> (a#(trans#0),a#(trans#1)))
    )

--permutePattern(patterns#10,perm3#3,perm3#4,perm2#0)

orbits = {};

--loop that applies all permutations to a zero pattern
--if the permuted zero pattern appears in the patterns list
--then it is removed and 1 is added to the orbit size
i = 0;
while i < #patterns do (
    orbitSize = 1;
    for rowPerm in perm3 do (
	for colPerm in perm3 do (
	    for trans in perm2 do (
	    	pattern = patterns#i;
	    	newPattern = permutePattern(pattern,rowPerm,colPerm,trans);
	    	newPattern = sort newPattern;
	    	if (newPattern != pattern and member(newPattern,patterns)) then (
		    patterns = delete(newPattern, patterns);
		    orbitSize = orbitSize + 1;
		    );
		);    
	    );
	);
    orbits = append(orbits,orbitSize);
    i = i + 1;
    )

#patterns
#orbits
sum orbits == #(subsets pos)

for i to #patterns - 1 do (
    print(i + 1);
    print(patterns#i);
    print(orbits#i);
    print("------------");
    )
