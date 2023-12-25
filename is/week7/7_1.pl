% サザエさんの人間関係
% 人(それぞれに性別を付与)
human(X) :- male(X).
human(X) :- female(X).
% 猫
cat(tama).
% 性別
female(fune).
female(wakame).
female(sazae).
female(hanazawa).
female(rika).
male(namihei).
male(katsuo).
male(nakajima).
male(masuo).
male(tarao).
% 配偶関係
wife(fune, namihei).
wife(sazae, masuo).
partner(X,Y) :- wife(X,Y). % XはYの配偶者
partner(X,Y) :- wife(Y,X). % XはYの配偶者
husband(X,Y) :- male(X),partner(X,Y) . % XはYの夫
% 親子関係
mother(fune, katsuo).
mother(fune, wakame).
mother(fune, sazae).
mother(sazae, tarao).
father(X,Y) :- husband(X,Z),mother(Z,Y) . % XはYの父
parent(X,Y) :- mother(X,Y). % XはYの親
parent(X,Y) :- father(X,Y). % XはYの親
child(X,Y) :- parent(Y,X). % XはYの子
% 友達関係
friend(katsuo,nakajima).
friend(katsuo,hanazawa).
friend(nakajima,hanazawa).
friend(tarao,rika).
friend(X,Y) :- friend(Y,X).
pet(tama,namihei).
% 一般論
siblings(X,Y) :- child(X,Z), child(Y,Z), X \= Y. % XとYは兄弟
uncle(X,Y) :- parent(Z,Y), siblings(X, Z), male(X). % XはYの叔父
aunt(X,Y) :- parent(Z,Y), siblings(X, Z), female(X). % XはYの叔母
ancestor(X,Y) :- parent(X,Y).
ancestor(X,Y) :- parent(X,Z), ancestor(Z,Y).
descendant(X,Y) :- ancestor(Y,X).

%ペットを一匹だけ飼っている人
one_pet(X) :- pet(Y,X), not((pet(Z,X), Z \= Y)).
