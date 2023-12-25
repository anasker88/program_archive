append([], Y, Y).
append([A|X], Y, [A|Z]) :- append(X, Y, Z).

pickoneup(A,[A|_]). %取り出した要素 、取り出し前のリスト
pickoneup(A,[_|Y]) :- pickoneup(A,Y).

remove(A,X,[A|X]).
remove(A,[B|X],[B|Y]) :- remove(A,X,Y).

%枝は[始点，終点]というリストのリスト。edgeは(始点、終点、枝のリスト)でそこに枝があるかを表す
edge(A,B,[[A,B]|_]). 
edge(A,B,[_|X]) :- edge(A,B,X).

%グラフ(V,E)上でSからTへのハミルトン路があるか。始点Sからある頂点Aへの枝があった時、VからSを除いた頂点集合Xについてグラフ(X,E)上でAからTへのハミルトン路を構成できることが条件。
makehamilton([S],_,S,S).
makehamilton(V,E,S,T) :- pickoneup(S,V),remove(S,X,V),edge(S,A,E),makehamilton(X,E,A,T).

hamilton(V,E) :- pickoneup(A,V),pickoneup(B,V),makehamilton(V,E,A,B).