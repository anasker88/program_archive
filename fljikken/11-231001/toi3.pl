% プレイヤーPの駒が三目並んでいるか
sanmoku(P,[[P,_,_],
           [P,_,_],
           [P,_,_]]).
sanmoku(P,[[_,P,_],
           [_,P,_],
           [_,P,_]]).
sanmoku(P,[[_,_,P],
           [_,_,P],
           [_,_,P]]).
sanmoku(P,[[P,P,P],
           [_,_,_],
           [_,_,_]]).
sanmoku(P,[[_,_,_],
           [P,P,P],
           [_,_,_]]).
sanmoku(P,[[_,_,_],
           [_,_,_],
           [P,P,P]]).
sanmoku(P,[[P,_,_],
           [_,P,_],
           [_,_,P]]).
sanmoku(P,[[_,_,P],
           [_,P,_],
           [P,_,_]]).
%一手先の行
next_line(P,[0,A,B],[P,A,B]).
next_line(P,[A,0,B],[A,P,B]).
next_line(P,[A,B,0],[A,B,P]).
%一手先の盤面
next_Banmen(P,[X,Y,Z],[W,Y,Z]) :- next_line(P,X,W).
next_Banmen(P,[X,Y,Z],[X,W,Z]) :- next_line(P,Y,W).
next_Banmen(P,[X,Y,Z],[X,Y,W]) :- next_line(P,Z,W).
%任意のプレイヤーについて一手先が存在しなければ終局
syuukyoku(B) :- \+ next_Banmen(_,B,_).
% 手番プレイヤーPが盤面Bで必勝であるか
win(P,_,B) :- sanmoku(P,B).%三目並んでいれば勝ち
win(_,Q,B) :- sanmoku(Q,B),!,false.%相手が三目並んでいればアウト
win(P,Q,B) :- next_Banmen(P,B,B_N),lose(Q,P,B_N).%一手打って相手が必ず負けなら必勝
%手番プレイヤーPが盤面Bで必ず負けか
lose(_,Q,B) :- sanmoku(Q,B).%相手が三目並んでいるなら負け
lose(P,Q,B) :- next_Banmen(P,B,B_N),  \+ win(Q,P,B_N),!,false.%ある手を打って相手が必勝でなくなるなら違う
lose(_,_,B) :- syuukyoku(B) , !, false.%そもそも終局なら違う
lose(_,_,_).%終局でなく、何を打っても相手が必勝なので負け
%必勝法が存在しない
no_winning_method :- \+ win(1,-1,[[0,0,0],[0,0,0],[0,0,0]]), \+ lose(1,-1,[[0,0,0],[0,0,0],[0,0,0]]).