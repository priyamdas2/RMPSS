
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]


NumericVector RMPSS_opt(NumericVector x0, Function func, double s_init = 1, double no_runs = 1000,
                        double max_iter = 10000, double rho_1 = 2, double rho_2 = 2, double phi = 10^(-20), double lambda = 10^(-10),
                        double tol_fun = 10^(-6), double tol_fun_2 = 10^(-20),int print_output = 0){
  double n = x0.size();
  NumericVector theta(n);
  double epsilon;
  double rho;
  NumericVector array_of_values(max_iter);
  NumericVector significant_positions(n);
  double num_of_significant;
  double epsilon_temp_pos;
  double epsilon_temp_neg;
  NumericVector possibility_pos(n);
  NumericVector possibility_neg(n);
  NumericVector total_lh_pos(n);
  NumericVector total_lh_neg(n);
  NumericMatrix matrix_update_at_h_pos(n,n);
  NumericMatrix matrix_update_at_h_neg(n,n);
  double current_lh;
  
  double garbage;
  
  
  double pos_movement_min;
  double neg_movement_min;
  int pos_movement_min_index;
  int neg_movement_min_index;
  
  double temp_sum;
  NumericVector new_theta(n);
  NumericVector yy;
  NumericMatrix theta_array(no_runs,n);
  
  
  // Initial point check if on simplex
  
  double sum_x0 = sum(x0);
  double min_x0 = min(x0);
  if(sum_x0 !=1 || min_x0 < 0)
  {Rcout << "ERROR! Provided starting point is not on simplex, its sum is : " << sum_x0 << std::endl;
    Rcout << "ERROR! Provided starting point is not on simplex, its min coordinate value is : " << min_x0 << std::endl;
    Rcout << "Default starting point is used, each coordinate is  :" << 1/n << std::endl;
    for(int dddd = 0; dddd < n; ++dddd)
    {x0[dddd] = 1/n;}}
  
  double start_value = Rcpp::as<double>(func(x0));
  
  for(int iii = 0; iii < no_runs; ++iii)
  {epsilon = s_init;
    if(iii == 0)
    {rho = rho_1;
      theta = clone(x0);
    }
    else
    {rho = rho_2;
      for(int gg = 0; gg < n; ++gg)
      {theta[gg] = theta_array(iii-1,gg);
      }
    }
    
    
    for(int i = 0; i < max_iter; ++i)
    {current_lh = Rcpp::as<double>(func(theta));
      
      if(print_output == 1)
      {Rcout << "Run number : " << iii+1 << std::endl;
        Rcout << "iteration : " << i+1 << std::endl;
        Rcout << "Current Obj function value :" << current_lh << std::endl;
        Rcout << "Current Theta sum :" << sum(theta) << std::endl;}
      
      
      
      // POSITIVE Movement
      
      for(int kkk = 0; kkk < n; ++kkk)
      {possibility_pos = clone(theta);
        
        garbage = 0;
        for(int lll = 0; lll < n; ++lll)
        {if(theta[lll] > lambda)
        {significant_positions[lll] = 1;
        }
        else
        {significant_positions[lll] = 0;
          garbage = garbage+theta[lll];}
        }
        significant_positions[kkk] = 0;
        
        if(theta[kkk] <= lambda)
        {garbage = garbage - theta[kkk];}
        
        
        num_of_significant = sum(significant_positions);
        
        
        
        if(num_of_significant == 0 || theta[kkk] == 1)
        {possibility_pos = clone(theta);}
        else
        {for(int lll = 0; lll < n; ++lll)
        {if(significant_positions[lll] == 1)
        {possibility_pos[lll] = theta[lll] - epsilon/num_of_significant + garbage/num_of_significant;
        }
        else
        {possibility_pos[lll] = 0;}
        }
        possibility_pos[kkk] = theta[kkk] + epsilon;
          
          
          epsilon_temp_pos = epsilon;
          
          while(min(possibility_pos) < 0 || max(possibility_pos)>1)
          {epsilon_temp_pos = epsilon_temp_pos/rho;
            possibility_pos[kkk] = theta[kkk] + epsilon_temp_pos;
            for(int lll = 0; lll < n; ++lll)
            {if(significant_positions[lll] == 1)
            {possibility_pos[lll] = theta[lll] - epsilon_temp_pos/num_of_significant + garbage/num_of_significant;
            }
            }
            
            if(epsilon_temp_pos < phi)
            {break;
            }
          }
        }
        
        if(min(possibility_pos)>= 0 && max(possibility_pos) <= 1)
        {total_lh_pos[kkk] = Rcpp::as<double>(func(possibility_pos));
        }
        else
        {possibility_pos = clone(theta);
          total_lh_pos[kkk] = current_lh;
        }
        
        
        
        for(int gg = 0; gg < n; ++gg)
        {matrix_update_at_h_pos(kkk,gg) = possibility_pos[gg];
        }
      }
      
      //NEGATIVE Movement
      
      for(int kkk = 0; kkk < n; ++kkk)
      {possibility_neg = clone(theta);
        
        
        garbage = 0;
        for(int lll = 0; lll < n; ++lll)
        {if(theta[lll] > lambda)
        {significant_positions[lll] = 1;
        }
        else
        {significant_positions[lll] = 0;
          garbage = garbage + theta[lll];}
        }
        significant_positions[kkk] = 0;
        
        if(theta[kkk] <= lambda)
        {garbage = garbage - theta[kkk];}
        
        num_of_significant = sum(significant_positions);
        if(num_of_significant == 0 || theta[kkk] == 1)
        {possibility_neg = clone(theta);}
        else
        {for(int lll = 0; lll < n; ++lll)
        {if(significant_positions[lll] == 1)
        {possibility_neg[lll] = theta[lll] + epsilon/num_of_significant + garbage/num_of_significant;
        }
        else
        {possibility_neg[lll] = 0;}
        }
        possibility_neg[kkk] = theta[kkk] - epsilon;
          
          
          epsilon_temp_neg = epsilon;
          
          while(max(possibility_neg) > 1 || possibility_neg[kkk] < 0)
          {epsilon_temp_neg = epsilon_temp_neg/rho;
            possibility_neg[kkk] = theta[kkk] - epsilon_temp_neg;
            for(int lll = 0; lll < n; ++lll)
            {if(significant_positions[lll] == 1)
            {possibility_neg[lll] = theta[lll] + epsilon_temp_neg/num_of_significant + garbage/num_of_significant;
            }
            }
            if(epsilon_temp_neg < phi)
            {break;
            }
          }
        }
        if(max(possibility_neg)<= 1 && possibility_neg[kkk]>=0)
        {total_lh_neg[kkk] = Rcpp::as<double>(func(possibility_neg));
        }
        else
        {possibility_neg = clone(theta);
          total_lh_neg[kkk] = current_lh;
        }
        
        for(int gg = 0; gg < n; ++gg)
        {matrix_update_at_h_neg(kkk,gg) = possibility_neg[gg];
        }
      }
      
      // Finding minimum index
      
      
      for(int gg = 0; gg < n; ++gg)
      {if(gg == 1)
      {pos_movement_min = total_lh_pos[gg];
        neg_movement_min = total_lh_neg[gg];
        pos_movement_min_index = 1;
        neg_movement_min_index = 1;
      }
      else
      {if(total_lh_pos[gg] < pos_movement_min)
      {pos_movement_min = total_lh_pos[gg];
        pos_movement_min_index = gg;
      }
      if(total_lh_neg[gg] < neg_movement_min)
      {neg_movement_min = total_lh_neg[gg];
        neg_movement_min_index = gg;
      }
      }
      }
      
      new_theta = theta;
      if(neg_movement_min < pos_movement_min)
      {if(neg_movement_min < current_lh)
        for(int kkk = 0; kkk<n;++kkk)
        {new_theta[kkk] = matrix_update_at_h_neg(neg_movement_min_index,kkk);
        }
      }
      else
      {if(pos_movement_min < current_lh)
        for(int kkk = 0; kkk<n;++kkk)
        {new_theta[kkk] = matrix_update_at_h_pos(pos_movement_min_index,kkk);
        }
      }
      
      theta = new_theta;
      
      array_of_values[i] = Rcpp::as<double>(func(theta));
      
      
      
      
      
      
      // Checking step-size
      
      if(i > 0)
      {if(abs(array_of_values[i] - array_of_values[i-1]) < tol_fun)
      {if(epsilon > phi)
      {epsilon = epsilon/rho;}
      else
      {break;
      }
      }
      }
      
    }  // end iterations of 1 run
    
    for(int kkk = 0; kkk<n;++kkk)
    {theta_array(iii,kkk) = theta[kkk];
    }
    
    if(iii > 0)
    {temp_sum = 0;
      for(int gg = 0; gg < n; ++gg)
      {temp_sum = temp_sum + pow(theta_array(iii,gg) - theta_array(iii-1,gg),2);
      }
      if(temp_sum/n < pow(tol_fun_2,2))
      {break;
      }
    }
    
  }
  Rcout << "Initial Obj function value :" << start_value << std::endl;
  Rcout << "Final Obj function value :" << current_lh << std::endl;
  return theta;
}