//Backend for the instanton finder in the AdS/Flat spacetime case.
#define DLL_EXPORT
#define INSTANTON_SOLVER_DLL_API //Comment this out if we wish to import these
//from a DLL instead.
#include "project_specific.h"

#include "multi_precision_definitions.h"
#include <boost/array.hpp>//For boost arrays (needed for ode_types)
#include "Flat.h"

//------------------------------------------------------------------------------
//Main function to export:
template< class value_type , class time_type , class solution_type,
          class solution_type_action >
void odeSolveFlatFixedBackground
    (value_type false_vacuum , value_type true_vacuum , value_type barrier ,
     potential< value_type >& V , time_type chimax , int odeSolverToUse,
     value_type RelTol, value_type AbsTol, value_type stepError, value_type xi,
     value_type lowerBound , value_type upperBound ,
     solution_grid< time_type , solution_type_action , value_type >& solOut,
     value_type& DSout,value_type precision,std::ostream& outStream,
     bool useSimpleBC,time_type epsilon_step)
{
//INPUTS:
    /*prhs[0] = a in GeV
     *prhs[1] = b (should be -1 < b < 1 to have a barrier)
     *prhs[2] = g (dimensionless coupling strength)
     *prhs[3] = chimax ("infinity" - ie, maximum time we allow)
     *prhs[4] = odeSolver, ie, which ode integration routine to use.
     *  odeSolver = 0 -> (Default) Cash-Karp method (Runge-Kutta)
     *  odeSolver = 1 -> Cash-Karp method (Runge-Kutta)
     *  odeSolver = 2 -> Dormand-Prince algorithm (Runge-Kutta)
     *  odeSolver = 3 -> Fehlberg algorithm (Runge-Kutta)
     *  odeSolver = 4 -> Bulirsch-Stoer algorithm (multi-step method)
     *prhs[5] = RelTol, optional (defaults to 1e-15)
     *prhs[6] = AbsTol, optional (defaults to 1e-15)
     *prhs[7] - stepError - optional. Fraction of Taylor series validity bound
     *          to use for first step. Defaults to 1e-15
     *prhs[8] - xi, non-minimal coupling. (optional - default to 0).
        NB - this won't actually affect the fixed background case in flat
            space!!!
     *prhs[9] - Lower bound on range to search (optional, defaults to barrier)
     *prhs[10] - Upper bound on range to search (optional, defaults to
                                                    true-vacuum)
     *prhs[11] - solOut - solution grid to store final result in.
     *prhs[12] - DSout - multi in which to store computed decay exponent.
     prhs[13] - precision - difference desired between upper and lower bounds
                on initial value of solution.
     */

    //Ok. we proceed by bisecting until we have a solution which overshoots
    //beyond
    //chimax, which stands in for infinity.

    //Overshoots occur when the solution crosses zero before chimax.
    //Undershoots when
    //its derivative crosses zero before then.

    //Instantiate the events function, state-type, dynamic BCs etc...
    //typedef boost::array< value_type , 2 > FlatSolType;//solution without
    //action, for speed.
    //typedef boost::array< value_type , 5 > AdSFlatSolTypeAction;
    typedef solution_type FlatSolType;
    typedef solution_type_action FlatSolTypeAction;
    typedef odeFlat< FlatSolType , time_type , value_type > odeRHSFlat;
    typedef eventsFlat< FlatSolType , time_type , value_type > evFlat;
    typedef dynamicBCsFlatFixed_taylor< value_type , FlatSolType , time_type >
            dynBCFlat;
    typedef odeFlatAction< FlatSolTypeAction , time_type , value_type >
            odeRHSFlatAction;
    typedef eventsFlat< FlatSolTypeAction , time_type , value_type >
            evFlatAction;
    typedef dynamicBCsFlatFixed_taylor< value_type , FlatSolTypeAction ,
                                        time_type > dynBCFlatAction;

    //Instantiation for overshoot/undershoot:
	value_type h = value_type(1.0);//h = phi_scale/M_p so h = 1 is Planck units.
	value_type V0 = value_type(0.0);//V at phi = 0.
    odeRHSFlat odeToSolve(V,xi, h, V0); //ode to solve (gravitational instanton
    //equation)
    //Events function - choose upper and lower bounds 5% outside false and true
    // vacuum, to leave us some leeway for numerical
    //fluctuations near to the nominal boundaries.
    evFlat events(false_vacuum,true_vacuum,false_vacuum,barrier,true_vacuum,
                  outStream); //events function.
        //Checks for events labelled IE:
        //IE = 1 -> y crosses zero (terminal event, => overshoot)
        //IE = 2 -> y' crosses zero (terminal event, => undershoot)
        //IE = 3 -> a crosses zero (non-terminal event, for reference only)
        //IE = 4 -> y crosses lowerBound (usually 0, but could be different -
                                          //terminal event)
        //IE = 5 -> y crosses upperBound (terminal event, indicates bad
                                          //initial condition or perhaps a->0)
    dynBCFlat bcAdjuster(V,stepError,1);//Performs analytic first step using a
                                        // taylor series approximation.
        //stepError - determines size of initial step.
        //version (last argument) - specifies which version of the ode we are
            //using.
            //version = 0 -> normal ode, no action.
            //version = 1 -> ode with action
            //version = 2 -> linearised ode.
    //bcFlat bcAdjusterSimple(V,xi,h,V0);



    //Search loop:
    int nMax = 200;
    int counter = 0;
    //Current bounds on location of bounce:
    value_type lower = lowerBound;//Barrier
    value_type upper = upperBound;//True vacuum
    //outStream << "Lower initial = " << lower << std::endl;
    //outStream << "Upper initial = " << upper << std::endl;
    value_type guess;
    multi diff = abs(upper - lower);
    while(diff > precision)
    {
        //Bisect to get initial guess:
        guess = (lower + upper)/multi(2.0);

        //Check for when the solution passes through its FWHM:
        events.value_to_check = (guess + false_vacuum)/value_type(2.0);

        //Formulate initial conditions vector:
        FlatSolType y0 = {{ guess , value_type(0.0) }};
        //Step away to avoid co-ordinate singularity at chi = 0:
        FlatSolType y01step;
        time_type epsilon;
        /*
        if(!useSimpleBC)
        {
            epsilon = bcAdjuster(y0,y01step,0);//Copies solution after step of
            //epsilon into y01step. Last argument tells us whether
        }
        else
        {
            epsilon = epsilon_step;
            y01step = bcAdjusterSimple(epsilon,y0);
        }
        */
        epsilon = bcAdjuster(y0,y01step);//Copies solution after step of
        //epsilon into y01step. Last argument tells us whether
            //we step forwards (0) or backwards (1) in time.
        //Set up time range:
        time_type tspanArray[2] = { epsilon , chimax };
        std::vector< time_type > tspan;
        tspan.assign(tspanArray,tspanArray + 2);
        //Solve the ode:
        solution_grid< time_type , FlatSolType , value_type > solTemp;
        time_type initStepMax = 0.01/time_type(chimax);
        int nSuccess = ode_solve<time_type,solution_type,value_type>
        (odeToSolve,y01step,tspan,events,odeSolverToUse,solTemp,RelTol,AbsTol,
         initStepMax);
        if(nSuccess != 1)
        {
            outStream << "nSuccess = " << nSuccess << " when solving for "
                      << "AdS-Flat instantons." << std::endl;
            throw "Integration failure.";
        }
        //Now check whether we found an overshoot or an undershoot:
        bool overshoots = false;
        bool undershoots = false;
        for(int j = 0; j < int(solTemp.IE.size());j++)
        {
            switch(solTemp.IE[j])
            {
            case 1:
                //overshoot.
                //outStream << "Overshoot." << std::endl;
                //reportToCaller("Overshoot.\n");
                overshoots = true;
                upper = guess;
                break;
            case 2:
                //undershoot, IF it occurs on the opposite side of the barrier
                //to the true vacuum (theoretically impossible to occur on the
                                      //other
                //side, so must be a numerical artefact if it does).
                if((true_vacuum - barrier)*(solTemp.YE[j][0] - barrier)
                   < value_type(0.0))
                {
                    //outStream << "Undershoot." << std::endl;
                    //reportToCaller("Undershoot.\n");
                    undershoots = true;
                    lower = guess;
                }
                break;
            case 3:
                //FWHM. No need to do anything.
                break;
            case 4:
                //Crossed lower bound. Indicates a bad initial condition.
                outStream << "Upper bound crossed during ode integration. "
                          << "Possibly a bad initial condition (wrong side of "
                          << "true vacuum) or solution undershot and event "
                          << "detection failed to detect this." << std::endl;
                //reportToCaller("Upper bound crossed during ode integration.
                //Possibly a bad initial condition (wrong side of true vacuum)
                //or solution undershot and event detection failed to detect
                //this.\n");
                throw "Unexpected ode-event encountered during integration.";
                break;
            case 5:
                //Crossed upper bound. Indicates a bad initial condition.
                outStream << "Lower bound crossed during ode integration. "
                          << "Possibly a bad initial condition (wrong side of "
                          << "true vacuum) or solution overshot and event"
                          << " detection failed to detect this." << std::endl;
                //reportToCaller("Upper bound crossed during ode integration.
                                 //Possibly a bad initial condition
                                 //(wrong side of true vacuum) or solution
                                 //undershot and event detection failed to
                                 //detect this.\n");
                throw "Unexpected ode-event encountered during integration.";
            default:
                outStream << "Unrecognised ode-event." << std::endl;
                //reportToCaller("Unrecognised ode-event.\n");
                throw "Unexpected ode-event encountered during integration.";
            }
            if(overshoots || undershoots)
            {
                break;
            }
        }
        counter++;
        diff = abs(lower - upper);


        if(overshoots)
        {
            //outStream << "Overshoot.\n" << std::endl;
        }
        else if(undershoots)
        {
            //outStream << "Undershoot.\n" << std::endl;
        }

        //Dump all available data if we get neither an overshoot nor an
        //undershoot, as this is
        //theoretically impossible so we want to understand what went wrong:
        else if(!(overshoots||undershoots))
        {
            outStream << "Neither undershoot nor overshoot detected. "
                      << "(Perhaps integration range should be extended?) "
                      << "Relevant data:\n";
            outStream << "no. of events = " << int(solTemp.IE.size())
                      << std::endl;
            for(int i = 0; i < int(solTemp.IE.size());i++)
            {
                outStream << "IE[" << i << "] = " << solTemp.IE[i] << std::endl;
                outStream << "YE[" << i << "] = " << std::endl;
                for(int j = 0;j < 2;j++)
                {
                    outStream << solTemp.YE[i][j] << std::endl;
                }
                outStream << "TE[" << i << "] = " << solTemp.TE[i] << std::endl;
            }
            outStream << "y0 = \n";
            for(int i = 0;i < 2;i++)
            {
                outStream << y01step[i] << std::endl;
            }
            outStream << "yend = \n";
            for(int i = 0;i < 2;i++)
            {
                outStream << solTemp.Y[int(solTemp.Y.size() - 1)][i]
                          << std::endl;
            }
            outStream << "Tend = " << solTemp.T[int(solTemp.T.size()) - 1]
                      << std::endl;
            //Now exit, throwing an error:
            throw "Error - Neither overshoot nor undershoot detected.\n";
        }

    }
    //We found the bounce. Are we within precision?
    if(diff > precision)
    {
        outStream << "Loop exited before finding bounce." << std::endl;
        //reportToCaller("Loop exited before finding bounce.\n");
    }

    //Compute action:
    odeRHSFlatAction odeToSolveAction(V,xi,value_type(1.0),value_type(0.0));
    evFlatAction eventsAction(false_vacuum,true_vacuum,false_vacuum,barrier,
                              true_vacuum , outStream);
    dynBCFlatAction bcAdjusterAction(V,stepError,2);

    //Check for when the solution passes through its FWHM:
    eventsAction.value_to_check = (guess + false_vacuum)/value_type(2.0);
    //Disable termination so that we can see the full extent of the solution:
    //eventsAction.terminate_on = false;
    //eventsAction.dynamicSwitch[0] = false;
    //eventsAction.dynamicSwitch[1] = false;

    FlatSolTypeAction y0Action = {{ guess , value_type(0.0) ,
                                    value_type(0.0) }};
    FlatSolTypeAction y0Action1step;
    time_type epsAction;
    /*
    if(!useSimpleBC)
    {
        epsAction = bcAdjusterAction(y0Action,y0Action1step,0);
        //Copies solution after step of epsilon into y01step.
        //Last argument tells us whether
    }
    else
    {
        epsAction = epsilon_step;
        y0Action1step = bcAdjusterSimple(epsAction,y0Action);
    }
    */
    epsAction = bcAdjusterAction(y0Action,y0Action1step);
    time_type tspanArray[2] = { epsAction , chimax };
    std::vector< time_type > tspan;
    tspan.assign(tspanArray,tspanArray + 2);
    time_type initStepMax = 0.01/time_type(chimax);
    int nSuccess = ode_solve<time_type,solution_type_action,value_type>
    (odeToSolveAction,y0Action1step,tspan,eventsAction,odeSolverToUse,solOut,
     RelTol,AbsTol,initStepMax);
    if(nSuccess != 1)
    {
        outStream << "nSuccess = " << nSuccess << " when solving for AdS-Flat "
                  << "instantons." << std::endl;
        throw "Integration failure.";
    }
    DSout = solOut.Y[int(solOut.T.size()) - 1][2];
    //Done.
}
//------------------------------------------------------------------------------
//Export an ode solver:
//Main function to export:
template< class value_type , class time_type , class solution_type,
          class solution_type_action >
void odeSolveFlatFixedBackgroundSingle
    (value_type y0,value_type false_vacuum , value_type true_vacuum ,
     value_type barrier , potential< value_type >& V , time_type chimax ,
     int odeSolverToUse, value_type RelTol, value_type AbsTol,
     value_type stepError, value_type xi, value_type lowerBound ,
     value_type upperBound , solution_grid< time_type , solution_type_action ,
     value_type >& solOut,value_type& DSout,value_type precision,
     std::ostream& outStream)
{
//INPUTS:
    /*prhs[0] = a in GeV
     *prhs[1] = b (should be -1 < b < 1 to have a barrier)
     *prhs[2] = g (dimensionless coupling strength)
     *prhs[3] = chimax ("infinity" - ie, maximum time we allow)
     *prhs[4] = odeSolver, ie, which ode integration routine to use.
     *  odeSolver = 0 -> (Default) Cash-Karp method (Runge-Kutta)
     *  odeSolver = 1 -> Cash-Karp method (Runge-Kutta)
     *  odeSolver = 2 -> Dormand-Prince algorithm (Runge-Kutta)
     *  odeSolver = 3 -> Fehlberg algorithm (Runge-Kutta)
     *  odeSolver = 4 -> Bulirsch-Stoer algorithm (multi-step method)
     *prhs[5] = RelTol, optional (defaults to 1e-15)
     *prhs[6] = AbsTol, optional (defaults to 1e-15)
     *prhs[7] - stepError - optional. Fraction of Taylor series validity bound
     *          to use for first step. Defaults to 1e-15
     *prhs[8] - xi, non-minimal coupling. (optional - default to 0).
        NB - this won't actually affect the fixed background case in
        flat space!!!
     *prhs[9] - Lower bound on range to search (optional, defaults to barrier)
     *prhs[10] - Upper bound on range to search (optional, defaults to
                                                 true-vacuum)
     *prhs[11] - solOut - solution grid to store final result in.
     *prhs[12] - DSout - multi in which to store computed decay exponent.
     prhs[13] - precision - difference desired between upper and lower bounds
        on initial value of solution.
     */

    //Ok. we proceed by bisecting until we have a solution which overshoots
    //beyond
    //chimax, which stands in for infinity.

    //Overshoots occur when the solution crosses zero before chimax.
    // Undershoots when
    //its derivative crosses zero before then.

    //Instantiate the events function, state-type, dynamic BCs etc...
    //typedef boost::array< value_type , 2 > FlatSolType;//solution without
    //action, for speed.
    //typedef boost::array< value_type , 5 > AdSFlatSolTypeAction;
    typedef solution_type FlatSolType;
    typedef solution_type_action FlatSolTypeAction;
    typedef odeFlat< FlatSolType , time_type , value_type > odeRHSFlat;
    typedef eventsFlat< FlatSolType , time_type , value_type > evFlat;
    typedef dynamicBCsFlatFixed_taylor< value_type , FlatSolType , time_type >
            dynBCFlat;
    typedef odeFlatAction< FlatSolTypeAction , time_type , value_type >
            odeRHSFlatAction;
    typedef eventsFlat< FlatSolTypeAction , time_type , value_type >
            evFlatAction;
    typedef dynamicBCsFlatFixed_taylor< value_type , FlatSolTypeAction ,
                                         time_type > dynBCFlatAction;

    //Instantiation for overshoot/undershoot:
    outStream << "Setup ode..." << std::endl;
	value_type h = value_type(1.0);//h = phi_scale/M_p so h = 1 is Planck units.
	value_type V0 = value_type(0.0);//V at phi = 0.
    odeRHSFlat odeToSolve(V,xi, h, V0); //ode to solve (gravitational
                                                        //instanton equation)

    //Compute action:
    odeRHSFlatAction odeToSolveAction(V,xi,value_type(1.0),value_type(0.0));
    evFlatAction eventsAction(false_vacuum*value_type(1.05),true_vacuum*
                              value_type(1.05),false_vacuum,barrier,
                              true_vacuum , outStream);
    dynBCFlatAction bcAdjusterAction(V,stepError,2);

    outStream << "Setup initial conditions..." << std::endl;
    FlatSolTypeAction y0Action = {{ y0 , value_type(0.0) , value_type(0.0) }};
    outStream << "V'(y0) = " << V.d(y0) << std::endl;
    FlatSolTypeAction y0Action1step;
    time_type epsAction = bcAdjusterAction(y0Action,y0Action1step);
    outStream << "y0Action1step = \n" << y0Action1step[0] << "\n"
              << y0Action1step[1] << "\n" << y0Action1step[2] << std::endl;
    outStream << "Setup timespan..." << std::endl;
    time_type tspanArray[2] = { epsAction , chimax };
    std::vector< time_type > tspan;
    tspan.assign(tspanArray,tspanArray + 2);
    outStream << "Integrating..." << std::endl;
    time_type initStepMax = 0.01/time_type(chimax);
    int nSuccess = ode_solve<time_type,solution_type_action,value_type>
    (odeToSolveAction,y0Action1step,tspan,eventsAction,odeSolverToUse,solOut,
     RelTol,AbsTol,initStepMax);
    if(nSuccess != 1)
    {
        outStream << "nSuccess = " << nSuccess << " when solving for AdS-Flat "
                  << "instantons." << std::endl;
        throw "Integration failure.";
    }
    outStream << "Compute action..." << std::endl;
    DSout = solOut.Y[int(solOut.T.size()) - 1][2];
    //Done.
    outStream << "Done." << std::endl;
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//Explicit instantiation and export of this function:
template DLL_EXPORT void odeSolveFlatFixedBackground
    < multi , multi , boost::array< multi , 2 > , boost::array< multi , 3 > >
    (multi , multi , multi , potential< multi >& , multi , int, multi, multi,
     multi, multi, multi , multi ,
     solution_grid< multi , boost::array< multi , 3 > , multi >& ,
     multi&,multi,std::ostream&,bool,multi);
//------------------------------------------------------------------------------
template DLL_EXPORT void odeSolveFlatFixedBackgroundSingle
    < multi , multi , boost::array< multi , 2 > , boost::array< multi , 3 > >
    (multi,multi , multi , multi , potential< multi >& , multi , int,
     multi, multi, multi, multi,multi , multi ,
     solution_grid< multi , boost::array< multi , 3 > , multi >& ,
     multi&,multi,std::ostream&);
//------------------------------------------------------------------------------
