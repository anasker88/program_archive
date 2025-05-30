/* toi1.pl */
male(kobo).
male(koji).
male(iwao).
female(sanae).
female(mine).
female(miho).

parent(kobo, koji).
parent(kobo, sanae).
parent(miho, koji).
parent(miho, sanae).
parent(sanae, iwao).
parent(sanae, mine).

father(X, Y) :- parent(X, Y), male(Y).
mother(X, Y) :- parent(X, Y), female(Y).

grandparent(X,Z) :- parent(X,Y), parent(Y,Z).

ancestor(X,Z) :- parent(X,Z).
ancestor(X,Z) :- parent(X,Y), ancestor(Y,Z).

sibling(X,Z) :- parent(X,Y), parent(Z,Y).

bloodrelative(X,Y) :- ancestor(X,Y).
bloodrelative(X,Y) :- ancestor(Y,X).
bloodrelative(X,Y) :- ancestor(X,Z) ,ancestor(Y,Z).
