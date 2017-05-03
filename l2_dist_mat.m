function D = l2_dist_mat(Y,X)

X_l2 = single((-2)*[ X ;ones( 1 , size(X,2) ); -1/2 * sum( X.^2 )  ]);

Y_l2 = single([ Y ; -1/2 * sum( Y.^2 ) ; ones( 1 , size(Y,2) ) ]);


D = cast(X_l2'*Y_l2,'like',Y);
