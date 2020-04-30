library(Rcpp)

library(RcppArmadillo)

cppFunction(depends="RcppArmadillo",plugins="cpp11",code="
            arma::vec Fbar(arma::vec& X, arma::vec& y) {
            arma::mat M=repmat(X,1,y.n_elem);
            arma::mat Y=repmat(trans(y),X.n_elem,1);
            arma::vec N(y.n_elem);
            N.ones();
            double n=X.n_elem;
            N=n*N;
            return(sum(trans(M>Y),1)/N);
            }
            ") #quantile c.s.f.
