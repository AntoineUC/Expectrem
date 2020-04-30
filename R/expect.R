library(Rcpp)

library(RcppArmadillo)

cppFunction(depends="RcppArmadillo",plugins="cpp11",code="
            arma::vec expect(arma::vec& X, arma::vec& alphas) {
            arma::mat M=repmat(X,1,alphas.n_elem);
            arma::vec a(arma::size(alphas));
            arma::vec b(arma::size(alphas));
            a.zeros();
            b.zeros();
            arma::vec gap(arma::size(alphas));
            gap.ones();
            arma::mat A=repmat(trans(a),X.n_elem,1);
            while(gap.max()>0.00001){
              A=repmat(trans(a),X.n_elem,1);
              b=mean(trans(M),1)-alphas%mean(trans(M-A)%trans(M<A),1)-(1-alphas)%mean(trans(M-A)%trans(M>A),1);
              gap=abs(b-a);
              a=b;
            }
            return(b);
            }
            ") #empirical expectile function
