#ifndef TRANSITION_SOLVER_H
#define TRANSITION_SOLVER_H
#include "ode_solver_generic_code.h"
//Attempts to solve the de-Sitter bounce problem by transitioning to an analytic
//form for a(x) (fixed background approximation) when entering the region where
//V0 dominates the scale factor.
//This should allow a more accurate computation of the action.

//The main workhorse of the solver should accept commands like:
/*
int nSuccess = ode_solve<time_type,dSSolTypeAction,value_type>
            (*odeToSolveAction,y0Action1step,tspan,*eventsAction,odeSolverToUse,
             LHS_solTemp,action_reltol,action_abstol,initStepMax);
*/
//This should hide the details of the switching between analytic regimes from
//dS.cpp. For all intents and purposes, it should look the same.
template<class time_type,class solution_type,class value_type>
int ode_solve_trans(odeRHS< solution_type, time_type, value_type >& odefun,
              solution_type& y0, std::vector< time_type >& tspan,
              events_function< solution_type, time_type, value_type >& events,
              int odeSolver,
              solution_grid< time_type, solution_type, value_type >& sol,
              value_type RelTol, value_type AbsTol,time_type initStepMaxFrac,
              value_type thresh,value_type y0_low,value_type y0_high,
              potential<value_type>& V,value_type V0,bool include_action,
              value_type xi,value_type h,bool true_vacuum_side,
              bool ap0termination,value_type fv)
{
    try
    {
    //Some needed constants
    const value_type _0p0 = value_type(0.0);
    const value_type _3p0 = value_type(3.0);
    const value_type _1p0 = value_type(1.0);
    const value_type _8p0 = value_type(8.0);
    const value_type _9p0 = value_type(9.0);
    const value_type _2p0 = value_type(2.0);
    const value_type _24p0 = value_type(24.0);
    const value_type PI = boost::math::constants::pi<value_type>();
    const value_type PI2 = PI*PI;
    const value_type H0 = h*sqrt(V0/(_3p0*(_1p0 - xi*fv*fv*h*h)));
    //Create a grid to store the temporary, fully-numerical solution:
    solution_grid< time_type , solution_type , value_type > tempGrid;

    //Setup new event detectors:
    //Version for the dynamic background method - terminates when |y'^2 + V|/V0
    // < thresh = RelTol so we can switch to the fixed background approximation:
    eventsSwitch< solution_type , time_type , value_type > evSwitch
        (y0_low,y0_high,V,thresh,V0);
    //Version for the fixed background approximation - terminates when
    // |y'^2 + V|/V0 > thresh so we can switch to the dynamic background
    //approximation.
    eventsSwitch_analytic< solution_type , time_type , value_type > evSwitchA
        (y0_low,y0_high,V,thresh,V0,h);
    //Setup analytic form for differential equation in the fixed background
    //approximation:
    ode_dS_switch< solution_type , time_type , value_type > odeToSolve
            (V,xi, h, V0);
    ode_dS_action_switch< solution_type , time_type , value_type > odeToSolveA
            (V,xi, h, V0,fv);

    //NB - ACTION CALCULATIONS ARE NOT YET SUITABLE FOR USE WITH NON-MINIMAL
    //COUPLING!!!!!



    //Turn off event termination for a' = 0 if requested. By default, the solver
    //terminates if we find a' = 0, which is desirable if we are computing from
    //both sides to compute the action, but undesirable if we are trying to
    //determine whether the solution is an overshoot or an undershoot:
    if(!ap0termination)
    {
        evSwitch.dynamicSwitch[0] = false;
        evSwitchA.dynamicSwitch[0] = false;
    }

    //Solve up to the first time we cross the threshold:
    int nSuccess;
    //Condition which decides whether to start by using the full numerical
    //solution, or the fixed background approximation:
    bool using_fixed_background = abs(V(y0[0]))/V0 < thresh;
    std::cout << "\nRatio = " << abs(V(y0[0]))/V0;
    if(include_action)
    {
        std::cout << "\nS_init = " << y0[4];
    }
    if(!using_fixed_background)
    {
        //Integrate with the dynamic background method. Boundary condition
        //setup should be done outside this function:
        nSuccess = ode_solve<time_type,solution_type,value_type>
                (odefun,y0,tspan,evSwitch,odeSolver,
                 tempGrid,RelTol,AbsTol,initStepMaxFrac);
    }
    else
    {
        //Integrate with the fixed background approximation. Have to setup the
        //boundary conditions, because the code outside this function assumes
        //that they are in the same form as above, which is not the case
        //when we want to treat a and a' analytically:
        std::vector<value_type> y0_fv_case;
        value_type S = include_action ? y0[4] : _0p0;
        value_type y0_fv_case_list[] = {y0[0],y0[1],S};
        if(include_action)
        {
            y0_fv_case.assign(y0_fv_case_list,y0_fv_case_list+3);
            nSuccess = ode_solve<time_type,solution_type,value_type>
                (odeToSolveA,y0_fv_case,tspan,evSwitchA,odeSolver,
                 tempGrid,RelTol,AbsTol,initStepMaxFrac);
        }
        else
        {
            y0_fv_case.assign(y0_fv_case_list,y0_fv_case_list+2);
            nSuccess = ode_solve<time_type,solution_type,value_type>
                (odeToSolve,y0_fv_case,tspan,evSwitchA,odeSolver,
                 tempGrid,RelTol,AbsTol,initStepMaxFrac);
        }
    }
    if(nSuccess != 1)
    {
        //Handle integration failures:
        std::cout << "nSuccess = " << nSuccess
                  << " when solving for AdS-Flat instantons."
                  << std::endl;
        sol.Y = tempGrid.Y;
        sol.T = tempGrid.T;
        sol.YE = tempGrid.YE;
        sol.TE = tempGrid.TE;
        sol.IE = tempGrid.IE;
        throw "Integration failure.";
    }

    //Check whether we ended by crossing the threshold. If not, then we are
    //done. If we crossed the threshold, then compute the second half, switching
    //from fixed background to dynamic background or vice versa:
    std::vector<value_type> varphi_list;//Phase shifts required to match a(t) at each
        //transition.
    std::vector<time_type> chi_0_list;//locations of the start of each fixed background
        //region.
    std::vector<time_type> chi_1_list;//locations of the start of each dynamic background
        //region.
    value_type phi = _0p0;//Initially, there isn't a phase shift.
    if(using_fixed_background)
    {
        chi_0_list.push_back(_0p0);
        chi_1_list.push_back(tempGrid.TE[int(tempGrid.IE.size()) - 1]);
        varphi_list.push_back(phi);
    }

    //Copy everything but the matching point from the old data
    // into the new solution (exclude the matching point to prevent it
    // occurring twice in the solution, since it is also the initial
    //condition for the second half - this is why we iterate only up to
    // int(tempGrid.Y.size()) - 1 rather than int(tempGrid.Y.size()) ):
    for(int i = 0;i < int(tempGrid.Y.size()) - 1;i++)
    {
        int nVecSize = include_action ? 5 : 4;
        std::vector<value_type> next_vector (nVecSize,_0p0);
        if(!using_fixed_background)
        {
            //Just copy the vector as is - it has all the expected
            //components already:
            next_vector = tempGrid.Y[i];
        }
        else
        {
            //Have to reconstruct a(t) using the analytic expression:
            value_type argi = H0*tempGrid.T[i] + phi;
            if(include_action)
            {
                value_type next_vec_list[] =
                    {tempGrid.Y[i][0],tempGrid.Y[i][1],
                     sin(argi)/H0,cos(argi),tempGrid.Y[i][2]};
                next_vector.assign(next_vec_list,next_vec_list + 5);
            }
            else
            {
                value_type next_vec_list[] =
                    {tempGrid.Y[i][0],tempGrid.Y[i][1],
                     sin(argi)/H0,cos(argi)};
                next_vector.assign(next_vec_list,next_vec_list + 4);
            }
        }
        /*
        if(i == 1)
        {
            std::cout << "\nnext_vector.size() = " << next_vector.size();
        }*/
        sol.Y.push_back(next_vector);
        /*
        if(int(sol.Y[i].size()) != 5)
        {
            std::cout << "\nError. Sol_Size = " << int(sol.Y[i].size());
        }*/
        sol.T.push_back(tempGrid.T[i]);
    }
    for(int i = 0;i < int(tempGrid.YE.size());i++)
    {
        //Same idea as the above loop, but for the events.
        int nVecSize = include_action ? 5 : 4;
        std::vector<value_type> next_vector (nVecSize,_0p0);
        if(!using_fixed_background)
        {
            next_vector = tempGrid.YE[i];
        }
        else
        {
            value_type argi = H0*tempGrid.T[i] + phi;
            if(include_action)
            {
                value_type next_vec_list[] =
                    {tempGrid.YE[i][0],tempGrid.YE[i][1],
                     sin(argi)/H0,cos(argi),tempGrid.YE[i][2]};
                next_vector.assign(next_vec_list,next_vec_list + 5);
            }
            else
            {
                value_type next_vec_list[] =
                    {tempGrid.YE[i][0],tempGrid.YE[i][1],
                     sin(argi)/H0,cos(argi)};
                next_vector.assign(next_vec_list,next_vec_list + 4);
            }
        }
        sol.YE.push_back(next_vector);
        sol.TE.push_back(tempGrid.TE[i]);
        sol.IE.push_back(tempGrid.IE[i]);
    }

    while(tempGrid.IE[int(tempGrid.IE.size()) - 1] == 6)
    {
        //Get the end of this solution as the new starting point:
        solution_type y_match = tempGrid.YE[int(tempGrid.IE.size()) - 1];
        time_type t_thresh = tempGrid.TE[int(tempGrid.IE.size()) - 1];
        value_type arg = H0*t_thresh + phi;
        /*value_type B_shift = _3p0*pi*pi*(_8p0 - _9p0*cos(arg) + cos(_3p0*arg))
                                        /(_2p0*V0);*/

        //Different solution type depending on whether we want to compute the
        //action or not:
        solution_type y_restart_action (3,_0p0);
        solution_type y_restart_no_action (2,_0p0);
        solution_type y_restart_action_fv_side (5,_0p0);
        solution_type y_restart_no_action_fv_side (4,_0p0);
        solution_type y_restart;
        //Assign different initial conditions for the second part, depending
        //on whether we are talking about a solution integrated from the LHS
        //or the RHS.
        if(include_action)
        {
            if(!using_fixed_background)
            {
                value_type y0_List_action[]
                = {y_match[0],y_match[1],y_match[4]};
                //Switch to the fixed background approximation:
                y_restart = y_restart_action;
                y_restart.assign(y0_List_action,y0_List_action + 3);
            }
            else
            {
                //Switch to the dynamic background:
                value_type y0_List_action_fv_side[]
                = {y_match[0],y_match[1],sin(arg)/H0,cos(arg),y_match[2]};
                y_restart = y_restart_action_fv_side;
                y_restart.assign(y0_List_action_fv_side,
                                 y0_List_action_fv_side + 5);
            }
        }
        else
        {
            if(!using_fixed_background)
            {
                //Switch to the fixed background approximation:
                value_type y0_List_no_action[] = { y_match[0] , y_match[1]};
                y_restart = y_restart_no_action;
                y_restart.assign(y0_List_no_action,
                                 y0_List_no_action + 2);
            }
            else
            {
                //Switch to the dynamic background method:
                value_type y0_List_no_action_fv_side[] =
                {y_match[0] , y_match[1], sin(arg)/H0, cos(arg)};
                y_restart = y_restart_no_action_fv_side;
                y_restart.assign(y0_List_no_action_fv_side,
                                 y0_List_no_action_fv_side + 4);
            }
        }

        //Now define a new grid to store the second half of the solution:
        solution_grid< time_type , solution_type , value_type > tempGrid2;


        if(include_action)
        {
            std::cout << "\nS_trans = " << y_restart[2];
        }
        for(int i = 0;i < int(y_restart.size());i++)
        {
            std::cout << "\ny_restart[" << i << "] = " << y_restart[i];
        }


        //Conduct the matching procedure to ensure the analytic
        //approximation matches the numerical solution at
        //the matching point:
        value_type a_end = using_fixed_background ?
                            sin(H0*t_thresh)/H0 : y_match[2];
        //std::cout << "\nH0*a_end = " << H0*a_end;
        //std::cout << "\nH0 = " << H0;
        //std::cout << "\na_end = " << a_end;
        phi = -H0*t_thresh + asin(H0*a_end);
        //std::cout << "\nphi = " << phi;
        if(!using_fixed_background)
        {
            varphi_list.push_back(phi);
        }
        odeToSolve.phi = phi;//(odeToSolve deduces H0 as part of
                                    //its constructor)
        odeToSolveA.phi = phi;
        //New tspan:
        std::vector< time_type > tspan_new = tspan;
        tspan_new[0] = t_thresh;
        //Proceed to solve from here onwards, only this time using the original
        //events function:
        if(include_action)
        {
            if(!using_fixed_background)
            {
                //Integrate with the fixed background approximation:
                nSuccess = ode_solve<time_type,solution_type,value_type>
                        (odeToSolveA,y_restart,tspan_new,evSwitchA,odeSolver,
                         tempGrid2,RelTol,AbsTol,initStepMaxFrac);
            }
            else
            {
                //Integrate with the dynamic background method:
                nSuccess = ode_solve<time_type,solution_type,value_type>
                        (odefun,y_restart,tspan_new,evSwitch,odeSolver,
                         tempGrid2,RelTol,AbsTol,initStepMaxFrac);
            }
        }
        else
        {
            if(!using_fixed_background)
            {
                //Integrate with the fixed background approximation:
                nSuccess = ode_solve<time_type,solution_type,value_type>
                        (odeToSolve,y_restart,tspan_new,evSwitchA,odeSolver,
                         tempGrid2,RelTol,AbsTol,initStepMaxFrac);
            }
            else
            {
                //Integrate with the dynamic background method:
                nSuccess = ode_solve<time_type,solution_type,value_type>
                        (odefun,y_restart,tspan_new,evSwitch,odeSolver,
                         tempGrid2,RelTol,AbsTol,initStepMaxFrac);
            }
        }
        if(nSuccess != 1)
        {
            std::cout << "nSuccess = " << nSuccess
                      << " when solving for AdS-Flat instantons."
                      << std::endl;
            sol.Y = tempGrid2.Y;
            sol.T = tempGrid2.T;
            sol.YE = tempGrid2.YE;
            sol.TE = tempGrid2.TE;
            sol.IE = tempGrid2.IE;
            throw "Integration failure.";
        }
        //Should now have a proper solution. Fill the output grid with it, using
        //an analytic approximation for a(t). Loop is the same as above with
        //tempGrid1, but which side uses the fixed background approximation
        //and which does not is reversed since we are swapping between the two
        //at the transition:
        for(int i = 0;i < int(tempGrid2.Y.size()) - 1;i++)
        {
            int nVecSize = include_action ? 5 : 4;
            std::vector<value_type> next_vector (nVecSize,_0p0);
            if(!using_fixed_background)
            {
                value_type argi = H0*tempGrid2.T[i] + phi;
                if(include_action)
                {
                    value_type next_vec_list[] =
                        {tempGrid2.Y[i][0],tempGrid2.Y[i][1],
                         sin(argi)/H0,cos(argi),tempGrid2.Y[i][2]};
                    next_vector.assign(next_vec_list,next_vec_list + 5);
                }
                else
                {
                    value_type next_vec_list[] =
                        {tempGrid2.Y[i][0],tempGrid2.Y[i][1],
                         sin(argi)/H0,cos(argi)};
                    next_vector.assign(next_vec_list,next_vec_list + 4);
                }
            }
            else
            {
                next_vector = tempGrid2.Y[i];
            }
            //std::cout << " " << int(next_vector.size());
            sol.Y.push_back(next_vector);
            sol.T.push_back(tempGrid2.T[i]);
        }
        for(int i = 0;i < int(tempGrid2.YE.size());i++)
        {
            int nVecSize = include_action ? 5 : 4;
            std::vector<value_type> next_vector (nVecSize,_0p0);
            if(!using_fixed_background)
            {
                value_type argi = H0*tempGrid2.T[i] + phi;
                if(include_action)
                {
                    value_type next_vec_list[] =
                        {tempGrid2.YE[i][0],tempGrid2.YE[i][1],
                         sin(argi)/H0,cos(argi),tempGrid2.YE[i][2]};
                    next_vector.assign(next_vec_list,next_vec_list + 5);
                }
                else
                {
                    value_type next_vec_list[] =
                        {tempGrid2.YE[i][0],tempGrid2.YE[i][1],
                         sin(argi)/H0,cos(argi)};
                    next_vector.assign(next_vec_list,next_vec_list + 4);
                }
            }
            else
            {
                next_vector = tempGrid2.YE[i];
            }
            //std::cout << " " << int(next_vector.size());
            sol.YE.push_back(next_vector);
            sol.TE.push_back(tempGrid2.TE[i]);
            sol.IE.push_back(tempGrid2.IE[i]);
        }
        //Setup the next round, if it is needed:
        if(!using_fixed_background)
        {
            chi_0_list.push_back(t_thresh);
            chi_1_list.push_back(tempGrid2.TE[int(tempGrid2.TE.size()) - 1]);
        }
        tempGrid = tempGrid2;
        using_fixed_background = !using_fixed_background;
    }

    //By this point, the solution has been successfully computed. Now we just
    //need to compute the residual analytic factor in the action, which should
    //mostly cancel (we hope).





    //Don't include the false vacuum part, otherwise it will be double
    //counted when we combine results from two integrals meeting in the middle.
    //All we want to sum up here are the accrued corrections to the action.
    value_type S_correction = _0p0;
    value_type A, B;
    //Residual factors to analytically subtract, corresponding
    //to the regions where we solve for a(t) analytically:
    for(int i = 0;i < int(chi_0_list.size());i++)
    {
        A = H0*chi_0_list[i] + varphi_list[i];
        B = H0*chi_1_list[i] + varphi_list[i];
        std::cout << "\nA = " << A;
        std::cout << "\nB = " << B;
        S_correction -= _3p0*PI2*(_9p0*cos(A) - cos(_3p0*A)
                                   - _9p0*cos(B) + cos(_3p0*B) )/_2p0;
    }
    S_correction /= V0;//Divide by V0 in Planck units.
    std::cout << "S_correction = " << S_correction;
    //Save the action correction so we can access it later:
    sol.S_correction = S_correction;

    //Add this correction onto the action of the solution:
    //Or don't. How would we separate it in the case where we integrate from
    //two directions and meet in the middle? Just save the result.
    /*
    std::cout << "\nS_correction = " << S_correction;
    for(int i = 0;i < int(sol.T.size());i++)
    {
        sol.Y[i][4] += S_correction;
    }*/


    }
    catch(char* c)
    {
        std::cout << c;
    }
}

#endif //TRANSITION_SOLVER_H
