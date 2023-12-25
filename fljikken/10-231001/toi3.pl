append([], Y, Y).
append([A|X], Y, [A|Z]) :- append(X, Y, Z).

reverse([],[]).
reverse([A|X],Y) :- reverse(X,REV_X),append(REV_X,[A],Y). 

concat([],[]).
concat([A|X],Y) :- append(A,Z,Y),concat(X,Z).