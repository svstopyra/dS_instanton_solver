//Backend for the instanton finder in the AdS/Flat spacetime case.
#define DLL_EXPORT
#define INSTANTON_SOLVER_DLL_API //Comment this out if we wish to import these
//from a DLL instead.

#include "project_specific.h"

#include "multi_precision_definitions.h"
#include <boost/array.hpp>//For boost arrays (needed for ode_types)
#include <boost/math/constants/constants.hpp>
#include "dS.h"
#include "find_zero.h"
#include "transition_solver.h"
#include <initializer_list>//Allow us to use initialiser lists for vector.
    //This simplifies the code somewhat.
#include "numerical_integration.h"
//------------------------------------------------------------------------------
//Struct to make passing constants easier:
template<class value_type, class time_type>
struct shooting_constants
{
public:
    const value_type stepError;
    const potential<value_type>& V;
    const time_type chimax;
    const bool linear_BCs;
    const bool use_analytic_a;
    const int odeSolverToUse;
    const time_type initStepMax;
    const bool reverse_shooting;


    //Constructor:
    shooting_constants(value_type STEPERROR,potential<value_type>& W,
                       time_type CHIMAX,bool LINEAR_BCS, bool USE_ANALYTIC_A,
                       int ODESOLVERTOUSE,time_type INITSTEPMAX,
                       bool REVERSE_SHOOTING)
    : stepError(STEPERROR), V(W), chimax(CHIMAX), linear_BCs(LINEAR_BCS),
      use_analytic_a(USE_ANALYTIC_A),odeSolverToUse(ODESOLVERTOUSE),
      initStepMax(INITSTEPMAX),reverse_shooting(REVERSE_SHOOTING)
    {
    }
};
//------------------------------------------------------------------------------
//Function which inspects a solution to determine whether it is an overshoot
//(true) or undershoot (false)
template< class value_type , class time_type >
bool check_overshoot
    (solution_grid< time_type , std::vector< value_type > ,
     value_type >& solTemp, value_type& upper,value_type& lower,
     const value_type& guess,std::ostream& outStream,int side,
     bool reverse_shoot,bool override_order = false);

//------------------------------------------------------------------------------
//Single shot, with initial conditions vector already given:
template< class time_type , class dSSolType , class value_type >
int doShot(value_type y0,dSSolType sol_initial,
           solution_grid<time_type,dSSolType,value_type>& solgrid,
           potential<value_type>& V,value_type V0,value_type xi,
           value_type h,value_type false_vacuum,
           shooting_constants<value_type,time_type>& constants,
           odeRHS<dSSolType,time_type,value_type>* odeToSolve,
           events_function< dSSolType , time_type, value_type >* events,
           value_type reltol,value_type abstol,value_type thresh,
           value_type overshoot_lower,value_type overshoot_upper,
           dynamicBCs< value_type , dSSolType , time_type >* bcAdjuster)
{
    //Clear grid:
    solgrid.erase();
    //Step away to avoid co-ordinate singularity at chi = 0:
    //outStream << "Loop stage 3" <<std::endl;
    dSSolType sol_initial_1step = sol_initial;
    time_type epsilon = (*bcAdjuster)(sol_initial,sol_initial_1step);
    for(int i = 0;i < sol_initial_1step.size();i++)
    {
        std::cout << "\nLHS_y01step[" << i << "] = " << sol_initial_1step[i];
    }

    //Time span:
    time_type tspanArray[2] = { epsilon , constants.chimax };
    std::vector< time_type > tspan;
    tspan.assign(tspanArray,tspanArray + 2);

    //Integration step:
    int nSuccess;
    if(constants.use_analytic_a)
    {
        nSuccess = ode_solve_trans<time_type,dSSolType,value_type>
              (*odeToSolve,sol_initial_1step,tspan,*events,
               constants.odeSolverToUse,solgrid,reltol,abstol,
               constants.initStepMax,thresh,overshoot_lower,overshoot_upper,V,
               V0,false,xi,h,false,false,false_vacuum);
    }
    else
    {
        nSuccess = ode_solve<time_type,dSSolType,value_type>
            (*odeToSolve,sol_initial_1step,tspan,*events,
             constants.odeSolverToUse,solgrid,reltol,abstol,
             constants.initStepMax);
    }
    if(nSuccess != 1)
    {
        std::cout << "nSuccess = " << nSuccess
                  << " when solving for AdS-Flat instantons."
                  << std::endl;
        throw "Integration failure.";
    }

    return nSuccess;
}

//------------------------------------------------------------------------------
//Single shot:
template< class time_type , class dSSolType , class value_type >
int doShotNoAction(value_type y0,
                   solution_grid<time_type,dSSolType,value_type>& solgrid,
                   potential<value_type>& V,value_type V0,value_type xi,
                   value_type h,value_type false_vacuum,
                   shooting_constants<value_type,time_type>& constants,
                   odeRHS<dSSolType,time_type,value_type>* odeToSolve,
                   events_function< dSSolType , time_type, value_type >* events,
                   value_type reltol,value_type abstol,value_type thresh,
                   value_type overshoot_lower,value_type overshoot_upper,
                   dynamicBCs< value_type , dSSolType , time_type >* bcAdjuster,
                   bool fixed_background)
{
    //Formulate initial conditions vector:
    dSSolType sol_initial;
    const value_type _0p0 = value_type(0.0);
    //value_type LHS_y0_List[] = { LHS_guess , _0p0 , _0p0 , a_scale };
    if(!fixed_background)
    {
        value_type sol_initial_List[] = { y0 , _0p0 , _0p0 , _0p0 };
        sol_initial.assign(sol_initial_List,sol_initial_List + 4);
    }
    else
    {
        value_type sol_initial_List[] = { y0 , _0p0 };
        sol_initial.assign(sol_initial_List,sol_initial_List + 2);
    }

    int nSuccess = doShot(y0,sol_initial,solgrid,V,V0,xi,h,false_vacuum,
                          constants,odeToSolve,events,reltol,abstol,thresh,
                          overshoot_lower,overshoot_upper,bcAdjuster);
    return nSuccess;
}
//------------------------------------------------------------------------------
//Single shot:
template< class time_type , class dSSolType , class value_type >
int doShotWithAction(value_type y0,
                   solution_grid<time_type,dSSolType,value_type>& solgrid,
                   potential<value_type>& V,value_type V0,value_type xi,
                   value_type h,value_type false_vacuum,
                   shooting_constants<value_type,time_type>& constants,
                   odeRHS<dSSolType,time_type,value_type>* odeToSolve,
                   events_function< dSSolType , time_type, value_type >* events,
                   value_type reltol,value_type abstol,value_type thresh,
                   value_type overshoot_lower,value_type overshoot_upper,
                   dynamicBCs< value_type , dSSolType , time_type >* bcAdjuster,
                   bool fixed_background)
{
    //Formulate initial conditions vector:
    dSSolType sol_initial;
    const value_type _0p0 = value_type(0.0);
    //value_type LHS_y0_List[] = { LHS_guess , _0p0 , _0p0 , a_scale };
    if(!fixed_background)
    {
        value_type sol_initial_List[] = { y0 , _0p0 , _0p0 , _0p0 ,_0p0, _0p0};
        sol_initial.assign(sol_initial_List,sol_initial_List + 6);
    }
    else
    {
        value_type sol_initial_List[] = { y0 , _0p0 , _0p0};
        sol_initial.assign(sol_initial_List,sol_initial_List + 3);
    }
    std::cout << "\nYes I am listening. This vector has " << sol_initial.size()
              << " components.";


    int nSuccess = doShot(y0,sol_initial,solgrid,V,V0,xi,h,false_vacuum,
                          constants,odeToSolve,events,reltol,abstol,thresh,
                          overshoot_lower,overshoot_upper,bcAdjuster);
    return nSuccess;
}
//------------------------------------------------------------------------------
//Function to print a list of events:
template<class value_type,class time_type,class dSSolType>
void printEvents(solution_grid< time_type , dSSolType , value_type >& solgrid)
{
    std::cout << "\nIE = (";
    for(int i = 0;i < solgrid.IE.size();i++)
    {
        if(i == 0)
        {
            std::cout << solgrid.IE[i];
        }
        else
        {
            std::cout << "," << solgrid.IE[i];
        }
    }
    std::cout << ").\n";
}
//------------------------------------------------------------------------------
//Check consistency of overshoot/undershoot after increasing precision:
template< class time_type , class dSSolType , class value_type >
void checkConsistency(std::vector<value_type>& y0_list,
                     std::vector<int>& end_type,
                     value_type& undershoot_y0,value_type& overshoot_y0,
                     solution_grid<time_type,dSSolType,value_type>& solgrid,
                     potential<value_type>& V,value_type V0,value_type xi,
                     value_type h,value_type false_vacuum,
                     shooting_constants<value_type,time_type>& constants,
                     odeRHS<dSSolType,time_type,value_type>* odeToSolve,
                     events_function<dSSolType,time_type,value_type>* events,
                     value_type reltol,value_type abstol,value_type thresh,
                     value_type overshoot_lower,value_type overshoot_upper,
                     dynamicBCs<value_type,dSSolType,time_type>* bcAdjuster,
                     value_type& bound_upper,value_type& bound_lower,int side,
                     bool using_action_version,bool fixed_background)
{
    std::cout << "\nChecking consistency for new precision.";
    bool lastUndershootValid = false;
    bool lastOvershootValid = false;
    bool undershoot_bound_updated = false;
    bool overshoot_bound_updated = false;
    int nSuccess;

    //Setup defaults, in case there are no overshoots or no undershoots for some
    //reason:
    bool reverse_shooting = constants.reverse_shooting;
    if((!reverse_shooting && side == 0) || reverse_shooting && side == 1)
    {
        overshoot_y0 = bound_lower;
        undershoot_y0 = bound_upper;
    }
    else if((!reverse_shooting && side==1) || (reverse_shooting && side==0))
    {
        undershoot_y0 = bound_lower;
        overshoot_y0 = bound_upper;
    }
    else
    {
        std::cout << "\nUnrecognised combination of 'side' "
                  << "and 'reverse_shooting'";
        throw "Invalid side and/or reverse_shooting.";
    }

    for(int i = y0_list.size() - 1;i >= 0;i--)
    {
        if( (lastUndershootValid && end_type[i] == 2) ||
            (lastOvershootValid && end_type[i] != 2 ) )
        {
            //Don't test if we already found a valid overshoot
            // or undershoot and the next element is of the same type as the
            //valid one - we want the one closest to the bounce, which will
            //be the first one we encounter.
            continue;
        }
        if(!using_action_version)
        {
            nSuccess = doShotNoAction(y0_list[i],solgrid,V,V0,
                      xi,h,false_vacuum,constants,odeToSolve,events,
                      reltol,abstol,thresh,
                      overshoot_lower,overshoot_upper,bcAdjuster,
                      fixed_background);
        }
        else
        {
            nSuccess = doShotWithAction(y0_list[i],solgrid,V,V0,
                      xi,h,false_vacuum,constants,odeToSolve,events,
                      reltol,abstol,thresh,
                      overshoot_lower,overshoot_upper,bcAdjuster,
                      fixed_background);
        }
        int new_end = solgrid.IE[solgrid.IE.size() - 1];
        //Update bounds, over-riding the previous requirement for new bounds to
        //be in order if we found an error.
        bool relevant_bound_updated = new_end == 2 ?
                                        undershoot_bound_updated :
                                            overshoot_bound_updated;
        //Override if the new end condition does not equal the old one, AND
        //We have not already found one of the old bounds for the relevant type
        //to be incorrect(if we have, don't override - just use the usual
        //ordering to update the bounds, because the latest bound is already
        //at the new precision - we only want to override if the current bound
        //was established at the old precision).
        /*
        //Was attempting to improve efficiency by updating the bounds each time
        //we compute a shot, but this actually causes issues in the search
        //algorithm with bounds failing to update because the guess falls
        //outside the bounds (likely due to the fact that the process of going
        //back to check doesn't compute the bisections 'in order').
        bool check = check_overshoot(solgrid,bound_upper,bound_lower,y0_list[i],
                                     std::cout,side,constants.reverse_shooting,
                                     (end_type[i] != new_end)
                                     && !relevant_bound_updated );*/
        if(!lastUndershootValid && end_type[i] == 2)
        {
            lastUndershootValid = (end_type[i] == 2 && new_end == 2);
            undershoot_bound_updated = !lastUndershootValid;
            if(lastUndershootValid)
            {
                std::cout << "\nFound valid undershoot.";
                undershoot_y0 = y0_list[i];
            }
            else
            {
                std::cout << "\nLast undershoot not valid."
                          << " Checking further back.";
            }
        }
        if(!lastOvershootValid && end_type[i] != 2)
        {
            lastOvershootValid = (end_type[i] != 2 && new_end != 2);
            overshoot_bound_updated = !lastOvershootValid;
            if(lastOvershootValid)
            {
                std::cout << "\nFound valid overshoot.";
                overshoot_y0 = y0_list[i];
            }
            else
            {
                std::cout << "\nLast overshoot not valid."
                          << " Checking further back.";
            }
        }
        if(lastUndershootValid && lastOvershootValid)
        {
            //We are done. Exit:
            std::cout << "Found valid undershoot and overshoot "
                      << "at new precision.";
            std::cout << "\nundershoot_y0 = " << undershoot_y0;
            std::cout << "\novershoot_y0 = " << overshoot_y0;
            break;
        }
    }
    if(!(lastUndershootValid && lastOvershootValid) )
    {
        std::cout << "\nCould not guarantee overshoot/undershoot bounds at new"
                  << " precision. Using widest available bounds.";
        std::cout << "\nundershoot_y0 = " << undershoot_y0;
        std::cout << "\novershoot_y0 = " << overshoot_y0;
    }
    //Update bounds with the newly found overshoots and undershoots:
    //bool reverse_shooting = constants.reverse_shooting;
    if((!reverse_shooting && side == 0) || reverse_shooting && side == 1)
    {
        bound_lower = overshoot_y0;
        bound_upper = undershoot_y0;
    }
    else if((!reverse_shooting && side==1) || (reverse_shooting && side==0))
    {
        bound_lower = undershoot_y0;
        bound_upper = overshoot_y0;
    }
    else
    {
        std::cout << "\nUnrecognised combination of 'side' "
                  << "and 'reverse_shooting'";
        throw "Invalid side and/or reverse_shooting.";
    }
    //Detect if the bounds switch order for some reason. This can happen for
    //various reasons, none of them good.
    if(bound_lower > bound_upper)
    {
        std::cout << "\nError. Bounds switched. Erroneous side "
                  << "variable, reverse_shooting, or error in"
                  << " check_overshoot.";
        throw "Bounds Switched!";
    }
}
//------------------------------------------------------------------------------
//Main function to export:
template< class value_type , class time_type >
void odeSolve_dS(value_type false_vacuum , value_type true_vacuum,
                 value_type barrier , potential< value_type >& V ,
                 time_type chimax_suggestion , int odeSolverToUse,
                 value_type RelTol, value_type AbsTol, value_type stepError,
                 value_type xi,
                 solution_grid< time_type , std::vector< value_type > ,
                 value_type >& solOut,
                 solution_grid< time_type , std::vector< value_type > ,
                 value_type >& solOutLHS,
                 solution_grid< time_type , std::vector< value_type > ,
                 value_type >& solOutRHS,
                 value_type& DSout,
                 value_type precision,value_type h,value_type W0,
                 std::vector<value_type>& other_return_values,
                 std::ostream& outStream,bool track_scale_factor,
                 bool linear_BCs,value_type action_reltol,
                 value_type action_abstol,bool graded_precision,
                 value_type a_scale,value_type overshoot_lower,
                 value_type overshoot_upper,
                 bool userSuppliedBounds,
                 value_type LHS_bound_lower,
                 value_type LHS_bound_upper,value_type RHS_bound_lower,
                 value_type RHS_bound_upper,
                 bool reverse_shooting,bool use_analytic_a,
                 bool fixed_background)
{
//INPUTS:

    try{
    //Ok. we proceed by bisecting until we have a solution which overshoots
    //beyond
    //chimax, which stands in for infinity.

    //Overshoots occur when the solution crosses zero before chimax.
    //Undershoots when
    //its derivative crosses zero before then.

    //Instantiate the events function, state-type, dynamic BCs etc...
    typedef std::vector< value_type > dSSolType;//solution without action,
        //for speed.
    typedef std::vector< value_type > dSSolTypeAction;
    typedef ode_dS_delta_a< dSSolType , time_type , value_type > odeRHS_dS;
    typedef events_delta_a< dSSolType , time_type , value_type > ev_dS;
    typedef dynamicBCsAdSFlat_taylor_delta_a< value_type , dSSolType,time_type >
            dynBC_dS;
    typedef ode_dS_action_delta_a< dSSolTypeAction , time_type , value_type >
            odeRHS_dSAction;
    typedef events_delta_a< dSSolTypeAction , time_type , value_type >
            ev_dSAction;
    typedef dynamicBCsAdSFlat_taylor_delta_a< value_type , dSSolTypeAction ,
                                      time_type > dynBCAdSFlatAction;

    typedef eventsGravFriction<dSSolType,time_type,value_type> eventsGrav;
    typedef odeGravFriction< dSSolType , time_type , value_type > odeGrav;
    typedef dynamicBCsGravFriction_taylor< value_type , dSSolType ,
                                      time_type > dynBCGrav;
    typedef odeGravFrictionAction<dSSolType,time_type,value_type>
                                    odeGravAction;
    //Fixed background approximation ode functions:
    typedef ode_dS_fixed< dSSolType , time_type , value_type > odeRHS_dS_fixed;
    typedef ode_dS_fixed_action< dSSolType , time_type , value_type >
            odeRHS_dS_fixed_action;
    typedef events_dS_fixed< dSSolTypeAction , time_type , value_type >
            ev_dS_fixed;
    typedef dynamicBCs_dS_fixed< value_type , dSSolTypeAction ,
                                      time_type > dynBC_dS_fixed;




    odeRHS_dS odeToSolveTrack(V,xi, h, W0); //ode to solve (gravitational
                                                           //instanton equation)
    odeRHS_dS_fixed odeToSolveFixed(V,xi, h, W0);
    odeGrav odeToSolveNoTrack(V,xi, h, W0);
    odeRHS<std::vector<value_type>,time_type,value_type>* odeToSolve;

    const value_type _2p0 = value_type(2.0);
    const value_type _0p0 = value_type(0.0);
    const value_type _1p0 = value_type(1.0);
    const value_type _6p0 = value_type(6.0);
    const value_type _4p0 = value_type(4.0);
    const value_type _8p0 = value_type(8.0);
    const value_type _3p0 = value_type(3.0);
    const value_type _0p01 = value_type(0.01);
    const value_type _24p0 = value_type(24.0);
    const value_type _10p0 = value_type(10.0);
    const value_type PI = boost::math::constants::pi<value_type>();
    const value_type eps = std::numeric_limits<value_type>::epsilon();

    //Determine whether we have user supplied bounds, and if not, refine the
    //locations of the critical points and use those instead:
    value_type LHS_lower;
    value_type LHS_upper;
    value_type RHS_lower;
    value_type RHS_upper;

    //Output precision:
    std::cout.precision(50);

    //RHS_upper = true_vacuum;
    //Determine whether or not to use the user supplied bounds:
    if(userSuppliedBounds)
    {
        LHS_lower = LHS_bound_lower;
        LHS_upper = LHS_bound_upper;
        RHS_lower = RHS_bound_lower;
        RHS_upper = RHS_bound_upper;
    }
    else
    {
        LHS_lower = false_vacuum;
        LHS_upper = barrier;
        RHS_lower = barrier;
        RHS_upper = true_vacuum;
    }


    /*
    value_type true_vacuum = tv;
    value_type barrier = bar;
    value_type LHS_lower = false_vacuum;
    value_type LHS_upper = bar;
    value_type RHS_lower = bar;
    value_type RHS_upper = tv;
    */
    std::cout << "\nfv = " << false_vacuum;
    std::cout << "\ntv = " << true_vacuum;

    //Specify the bounds used to identify overshoots (defaults to the false
    // and true vacuum)
    value_type H0 = h*sqrt(W0/_3p0);
    ev_dS eventsTrack(overshoot_lower,overshoot_upper,H0);
    ev_dS_fixed eventsdSFixed(overshoot_lower,overshoot_upper,H0);
    //events function.
    eventsGrav eventsNoTrack(overshoot_lower,overshoot_upper);
    events_function< std::vector<value_type> , time_type, value_type >* events;
    //ev_dS events(false_vacuum,true_vacuum); //events function.
        //Checks for events labelled IE:
        //IE = 1 -> a' crosses zero (terminal event, => overshoot)
        //IE = 2 -> y' crosses zero (terminal event, => undershoot)
        //IE = 3 -> a crosses zero (non-terminal event, for reference only)
        //IE = 4 -> y crosses upperBound (terminal event, indicates overshoot)
        //IE = 5 -> y crosses lowerBound (terminal event, indicates a->0)

    //Note - switch off termination for the event where we detect a'(x) = 0.
    //This will hinder the overshoot/undershoot method if we leave it on
    //because sometimes it isn't possible to tell whether a solution is an
    //overshoot or undershoot just by looking at the solution up until
    //a' = 0. Leave it on for the action version, however, because we only need
    //to compute half the solution from either side in that case anyway:
    eventsTrack.dynamicSwitch[0] = false;//Switches off termination for event 1
    eventsNoTrack.dynamicSwitch[0] = false;
    eventsdSFixed.dynamicSwitch[0] = false;


    //Setup the functions which handle the boundary conditions at x = 0. This is
    //a co-ordinate singularity, so we need to perform a small step away from
    //it using a Taylor expansion in order to avoid numerical difficulties;

    //Two different versions, depending on whether we track a, a' separately,
    //('track' version) or just use an equation for a'/a ('no track' version).
    //Final argument determines whether we are using boundary conditions for the
    //action or not ('1' => no action, '2' => computing action)
    dynBC_dS bcAdjusterTrack(xi,h,W0,V,stepError,1);
    dynBCGrav bcAdjusterNoTrack(xi,h,W0,V,stepError,1);
    dynBC_dS_fixed bcAdjuster_ds_fixed(xi,h,W0,V,stepError,1);
    dynamicBCs_basicMethod< value_type, dSSolType, time_type >
        dynBCvac(xi,h,W0,V,stepError);//Applies to solutions computed without
            //action.
    dynamicBCs_actionMethod< value_type, dSSolTypeAction, time_type >
        dynBCvacAction(xi,h,W0,V,stepError);//applies to solutions computed
            //with action.
    //Setup pointers to the correct boundary conditions function, ode,
    //and events function, depending on whether
    //we use the 'track' or 'no track' version.
    dynamicBCs< value_type , std::vector<value_type> , time_type >* bcAdjuster;

    if(fixed_background)
    {
        std::cout << "\nUsing fixed background approximation: ";
        bcAdjuster = &bcAdjuster_ds_fixed;
        odeToSolve = &odeToSolveFixed;
        events = &eventsdSFixed;
    }
    else if(track_scale_factor)
    {
        std::cout << "\nTracking full scale factor.";
        bcAdjuster = &bcAdjusterTrack;
        odeToSolve = &odeToSolveTrack;
        events = &eventsTrack;
    }
    else
    {
        std::cout << "\nUsing a'/a only.";
        bcAdjuster = &bcAdjusterNoTrack;
        odeToSolve = &odeToSolveNoTrack;
        events = &eventsNoTrack;
    }
    //Performs analytic first
    //step using a taylor series approximation.
        //stepError - determines size of initial step.
        //version (last argument) - specifies which version of the ode we are
            //using.
            //version = 1 -> normal ode, no action.
            //version = 2 -> ode with action

    //Analogous functions for the action version:
    //outStream << "\nSetup action version.";
    odeRHS_dSAction odeToSolveActionTrack(V,xi,h,W0);
    odeRHS_dS_fixed_action odeToSolveFixedAction(V,xi, h, W0);
    //outStream << "\n1";
    ev_dSAction eventsActionTrack(false_vacuum,_2p0*true_vacuum,H0);
    eventsActionTrack.dynamicSwitch[0] = false;
    ev_dS_fixed eventsActionFixed(false_vacuum,_2p0*true_vacuum,H0);
    eventsActionFixed.dynamicSwitch[0] = false;
    //outStream << "\n2";
    dynBCAdSFlatAction bcAdjusterActionTrack(xi,h,W0,V,stepError,2);
    dynBC_dS_fixed bcAdjusterActionFixed(xi,h,W0,V,stepError,2);
    //outStream << "\n3";
    odeGravAction odeToSolveActionNoTrack(V,xi,h,W0);
    //outStream << "\n4";
    eventsGrav eventsActionNoTrack(false_vacuum,_2p0*true_vacuum);
    eventsActionNoTrack.dynamicSwitch[0] = false;
    //outStream << "\n5";
    dynBCGrav bcAdjusterActionNoTrack(xi,h,W0,V,stepError,2);
    //outStream << "\n6";

    dynamicBCs< value_type , std::vector<value_type> , time_type >*
        bcAdjusterAction;
    //outStream << "\n7";
    events_function< std::vector<value_type> , time_type, value_type >*
        eventsAction;
    //outStream << "\n8";
    odeRHS<std::vector<value_type>,time_type,value_type>* odeToSolveAction;
    //outStream << "\n9";
    if(fixed_background)
    {
        odeToSolveAction = &odeToSolveFixedAction;
        bcAdjusterAction = &bcAdjusterActionFixed;
        eventsAction = &eventsActionFixed;
    }
    else if(track_scale_factor)
    {
        odeToSolveAction = &odeToSolveActionTrack;
        bcAdjusterAction = &bcAdjusterActionTrack;
        eventsAction = &eventsActionTrack;
    }
    else
    {
        odeToSolveAction = &odeToSolveActionNoTrack;
        bcAdjusterAction = &bcAdjusterActionNoTrack;
        eventsAction = &eventsActionNoTrack;
    }
    //outStream << "\n10";

    //Setup bounds to avoid integrating in a region we won't actually need:
    /*eventsActionTrack.y0upperBound = RHS_upper;
    eventsActionTrack.y0lowerBound = LHS_lower;
    eventsActionNoTrack.y0upperBound = RHS_upper;
    eventsActionNoTrack.y0lowerBound = LHS_lower;*/
    //Keep termination for a' = 0 event on in this case - we only need half
    //the solution from either side!


    //Check for when the solution passes through its FWHM:
    //eventsAction.value_to_check = (guess + false_vacuum)/value_type(2.0);
    //std::cout << "\nSetup bounce BCs.";





    //Search loop:
    int nMax = 200;//Terminate after this many iterations and generate a
        //a warning if we haven't found the bounce.
    int counter = 0;




    //Output some diagnostic parameters:
    //outStream << "Lower initial = " << lower << std::endl;
    //outStream << "Upper initial = " << upper << std::endl;

    outStream << "\nbarrier = " << barrier;
    outStream << "\nfv = " << false_vacuum;
    outStream << "\ntv = " << true_vacuum;

    outStream << "\nLHS_upper = " << LHS_upper;
    outStream << "\nRHS_upper = " << RHS_upper;
    outStream << "\nLHS_lower = " << LHS_lower;
    outStream << "\nRHS_lower = " << RHS_lower;


    //Difference between the bounds on the bounce (gives an estimate of how
    //cllose we are):
    value_type LHS_diff = abs(LHS_upper - LHS_lower);
    value_type RHS_diff = abs(RHS_upper - RHS_lower);
    outStream << "LHS_diff = " << LHS_diff << std::endl;
    outStream << "RHS_diff = " << RHS_diff << std::endl;
    //Initial guess at bounce location:
    value_type LHS_guess = (LHS_lower + LHS_upper)/_2p0;
    value_type RHS_guess = (RHS_lower + RHS_upper)/_2p0;


    //Tolerances to decide whether we are sufficiently close to the bounce
    //that we can stop looking:
    //value_type tol_LHS = precision*(abs(LHS_guess) + _1p0)/_2p0;
    //value_type tol_RHS = precision*(abs(RHS_guess) + _1p0)/_2p0;
    value_type tol_LHS = precision*abs(LHS_guess);
    value_type tol_RHS = precision*abs(RHS_guess);
    outStream << "tol_LHS = " << tol_LHS << std::endl;
    outStream << "tol_RHS = " << tol_RHS << std::endl;

    //Maximum value of x to integrate up to:
    time_type chimax = chimax_suggestion;//Can add code here to adjust this
        //dynamically if we want, ie, take the user's input as a suggestion
        //only. May not be necessary, however.
    std::cout << "\nchimax = " << chimax;
    time_type initStepMax = time_type(0.01)/chimax;
    //Checks whether a(x) matches for the LHS and RHS solutions.
    bool a_match = false;

    //If we choose to steadily ramp up the precision as we get closer to the
    //bounce, use these initial values for the ode tolerances:
    value_type graded_reltol = RelTol;
    value_type graded_abstol = AbsTol;

    //structs to store the temporary results:
    solution_grid< time_type , dSSolType , value_type > LHS_solTemp;
    solution_grid< time_type , dSSolType , value_type > RHS_solTemp;

    //Storage for initial conditions:
    std::vector<value_type> y0_rhs_list;
    y0_rhs_list.push_back(RHS_upper);
    y0_rhs_list.push_back(RHS_lower);
    std::vector<value_type> y0_lhs_list;
    y0_lhs_list.push_back(LHS_upper);
    y0_lhs_list.push_back(LHS_lower);
    std::vector<int> end_type_rhs;
    end_type_rhs.push_back(reverse_shooting ? 2 : 5);
    end_type_rhs.push_back(reverse_shooting ? 5 : 2);
    std::vector<int> end_type_lhs;
    end_type_lhs.push_back(reverse_shooting ? 5 : 2);
    end_type_lhs.push_back(reverse_shooting ? 2 : 5);

    value_type thresh = precision;
    //For handling graded precision:
    value_type grading_thresh_lhs = value_type(1e-5);
    value_type grading_thresh_rhs = value_type(1e-5);
    value_type grade = value_type(1e-5);
    value_type graded_reltol_lhs = value_type(1e-10);
    value_type graded_reltol_rhs = graded_reltol_lhs;
    value_type graded_abstol_lhs = graded_abstol;
    value_type graded_abstol_rhs = graded_abstol;
    value_type last_undershoot_lhs;
    value_type last_overshoot_lhs;
    value_type last_undershoot_rhs;
    value_type last_overshoot_rhs;
    if(!reverse_shooting)
    {
        last_undershoot_lhs = LHS_upper;
        last_undershoot_rhs = RHS_lower;
        last_overshoot_lhs = LHS_lower;
        last_overshoot_rhs = RHS_upper;
    }
    else
    {
        last_undershoot_lhs = LHS_lower;
        last_undershoot_rhs = RHS_upper;
        last_overshoot_lhs = LHS_upper;
        last_overshoot_rhs = RHS_lower;
    }

    //Struct to make passing constants easier:
    shooting_constants<value_type,time_type> constants
    (stepError,V,chimax,linear_BCs, use_analytic_a,odeSolverToUse,initStepMax,
     reverse_shooting);

    int undershoots_lhs_found = 0;
    int overshoots_lhs_found = 0;
    int undershoots_rhs_found = 0;
    int overshoots_rhs_found = 0;


    //Bounce search loop:
    bool inconsistent = true;
    bool inconsistent_lhs = true;
    bool inconsistent_rhs = true;
    bool compute_action = false;
    bool both_undershoots = false;
    bool found_solution_lhs = false;
    bool found_solution_rhs = false;
    bool found_solution = false;
    bool no_undershoots_found_lhs = true;
    bool no_undershoots_found_rhs = true;
    do
    {
        //Integrate LHS solution (but don't bother if we already reached max
        //precision on it - just skip it and wait for the RHS to catch up).
        //solution_grid< time_type , dSSolType , value_type > LHS_solTemp;
        if(!found_solution_lhs)
        {
            LHS_solTemp.erase();
            //outStream << "Loop stage 1" <<std::endl;
            //Bisect to get initial guess:
            LHS_guess = (LHS_lower + LHS_upper)/_2p0;
            //tol_LHS = precision*(abs(LHS_guess) + _1p0)/_2p0;
            //tol_LHS = precision*abs(LHS_guess);
            tol_LHS = precision;
            outStream << "\nLHS_diff = " << LHS_diff;
            outStream << "\nLHS_y0 = " << LHS_guess;
            //outStream << "\nV'(LHS_y0) = " << V.d(LHS_guess);


            //Solve ode:
            outStream << "\nIntegrating LHS..." <<std::endl;
            int LHS_success;
            if(graded_precision)
            {
                if( (abs(LHS_diff) < grading_thresh_lhs) &&
                     (undershoots_lhs_found > 3 && overshoots_lhs_found > 3))
                {
                    //Old version of graded precision - suspect this doesn't
                    // work
                    //this well as it may mis-identify solutions:
                    /*
                    graded_reltol = std::min(value_type(1e-5)*
                                             abs(LHS_diff/LHS_guess),RelTol);
                    graded_abstol = AbsTol;
                    */
                    //New version - only increase the precision every so often,
                    // and
                    //when we do this, check the last two solutions are the same
                    //nature at the new precision:

                    //Upgrade precision:
                    grading_thresh_lhs *= grade;
                    graded_reltol_lhs *= grade;
                    graded_abstol_lhs = AbsTol;
                    //Check that the last overshoot/undershoot remains an
                    //overshoot/undershoot:
                    checkConsistency(y0_lhs_list,end_type_lhs,
                                     last_undershoot_lhs,last_overshoot_lhs,
                                     LHS_solTemp,V,W0,xi,h,false_vacuum,
                                     constants,
                                     odeToSolve,events,graded_reltol_lhs,
                                     graded_abstol_lhs,thresh,overshoot_lower,
                                     overshoot_upper,bcAdjuster,LHS_upper,
                                     LHS_lower,0,compute_action,
                                     fixed_background);
                    //Modify the guess accordingly:
                    LHS_guess = (last_overshoot_lhs + last_undershoot_lhs)/_2p0;
                    std::cout << "\nLHS_upper = " << LHS_upper;
                    std::cout << "\nLHS_lower = " << LHS_lower;
                    std::cout << "\nLHS_diff = " << LHS_upper - LHS_lower;
                    /*
                    if(!reverse_shooting)
                    {
                        LHS_upper = last_undershoot_lhs;
                        LHS_lower = last_overshoot_lhs;
                        LHS_diff = abs(LHS_upper - LHS_lower);
                    }
                    else
                    {
                        LHS_upper = last_overshoot_lhs;
                        LHS_lower = last_undershoot_lhs;
                        LHS_diff = abs(LHS_upper - LHS_lower);
                    }*/
                }
            }
            else
            {
                if(inconsistent_lhs)
                {
                    graded_reltol_lhs = RelTol;
                    graded_abstol_lhs = AbsTol;
                }
                else
                {
                    graded_reltol_lhs = action_reltol;
                    graded_abstol_lhs = action_abstol;
                }
            }
            std::cout << "\ngraded_reltol = " << graded_reltol_lhs;
            //Perform shot:
            if(!compute_action)
            {
                //Version without action:
                LHS_success = doShotNoAction(LHS_guess,LHS_solTemp,V,W0,xi,
                                   h,false_vacuum,constants,odeToSolve,events,
                                   graded_reltol_lhs,graded_abstol_lhs,thresh,
                                   overshoot_lower,overshoot_upper,bcAdjuster,
                                   fixed_background);
            }
            else
            {
                //Version with action:
                LHS_success = doShotWithAction(LHS_guess,LHS_solTemp,V,W0,xi,
                                   h,false_vacuum,constants,odeToSolveAction,
                                   eventsAction,action_reltol,action_abstol,
                                   thresh,overshoot_lower,overshoot_upper,
                                   bcAdjusterAction,
                                   fixed_background);
            }


            //Now check whether we found an overshoot or an undershoot:
            bool LHS_overshoots = check_overshoot(LHS_solTemp,LHS_upper,
                                                  LHS_lower,
                                                  LHS_guess,outStream,0,
                                                  reverse_shooting);
            if(LHS_overshoots)
            {
                overshoots_lhs_found++;
            }
            else
            {
                undershoots_lhs_found++;
            }
            std::cout << "\n" << (LHS_overshoots ? "overshoot " : "undershoot ")
                      << "\nLHS_lower = " << LHS_lower
                      << "\nLHS_upper = " << LHS_upper
                      << "\nLHS_guess = " << LHS_guess;
            printEvents(LHS_solTemp);
            end_type_lhs.push_back(LHS_solTemp.IE[LHS_solTemp.IE.size() - 1]);
            y0_lhs_list.push_back(LHS_guess);
            //Record whether we have found the solution. Only accept it
            //if we have better precision than the tolerance and
            //have verified that the action has been computed:
            found_solution_lhs = (LHS_diff < tol_LHS) && compute_action
                                    && (end_type_lhs.back() == 2);
        }



        //Integrate RHS solution (but don't bother if we already reached max
        //precision on it - just skip it and wait for the LHS to catch up).
        //solution_grid< time_type , dSSolType , value_type > RHS_solTemp;
        if(!found_solution_rhs)
        {
            RHS_solTemp.erase();
            RHS_guess = (RHS_lower + RHS_upper)/_2p0;
            //tol_RHS = precision*(abs(RHS_guess) + _1p0)/_2p0;
            //tol_RHS = precision*abs(RHS_guess);
            tol_RHS = precision;
            outStream << "\nRHS_diff = " << RHS_diff;
            outStream << "\nRHS_y0 = " << RHS_guess;
            //outStream << "\nV'(RHS_y0) = " << V.d(RHS_guess);
            //outStream << "Loop stage 2" <<std::endl;
            //Search for FWHM
            //events.value_to_check = (RHS_guess + LHS_guess)/_2p0;

            //Solve the ode:
            outStream << "\nIntegrating RHS..." <<std::endl;
            if(graded_precision)
            {
                if( (abs(RHS_diff) < grading_thresh_rhs) &&
                 (undershoots_rhs_found > 3 && overshoots_rhs_found > 3) )
                {
                    //Old version of graded precision - suspect this doesn't
                    // work
                    //this well as it may mis-identify solutions:
                    /*
                    graded_reltol = std::min(value_type(1e-5)*
                                             abs(LHS_diff/LHS_guess),RelTol);
                    graded_abstol = AbsTol;
                    */
                    //New version - only increase the precision every so often,
                    // and
                    //when we do this, check the last two solutions are the same
                    //nature at the new precision:

                    //Upgrade precision:
                    grading_thresh_rhs *= grade;
                    graded_reltol_rhs *= grade;
                    graded_abstol_rhs = AbsTol;
                    //Check that the last overshoot/undershoot remains an
                    //overshoot/undershoot:
                    checkConsistency(y0_rhs_list,end_type_rhs,
                                     last_undershoot_rhs,last_overshoot_rhs,
                                     RHS_solTemp,V,W0,xi,h,false_vacuum,
                                     constants,
                                     odeToSolve,events,graded_reltol_rhs,
                                     graded_abstol_rhs,thresh,overshoot_lower,
                                     overshoot_upper,bcAdjuster,RHS_upper,
                                     RHS_lower,1,compute_action,
                                     fixed_background);
                    //Modify the guess accordingly:
                    RHS_guess = (last_overshoot_rhs + last_undershoot_rhs)/_2p0;
                    std::cout << "\nRHS_upper = " << RHS_upper;
                    std::cout << "\nRHS_lower = " << RHS_lower;
                    std::cout << "\nRHS_diff = " << RHS_upper - RHS_lower;
                    /*
                    if(!reverse_shooting)
                    {
                        RHS_upper = last_overshoot_rhs;
                        RHS_lower = last_undershoot_rhs;
                        RHS_diff = abs(RHS_upper - RHS_lower);
                    }
                    else
                    {
                        RHS_upper = last_undershoot_rhs;
                        RHS_lower = last_overshoot_rhs;
                        RHS_diff = abs(RHS_upper - RHS_lower);
                    }*/
                }
            }
            else
            {
                if(inconsistent_rhs)
                {
                    graded_reltol_rhs = RelTol;
                    graded_abstol_rhs = AbsTol;
                }
                else
                {
                    graded_reltol_rhs = action_reltol;
                    graded_abstol_rhs = action_abstol;
                }
            }
            std::cout << "\ngraded_reltol = " << graded_reltol_rhs;
            int RHS_nSuccess;
            if(!compute_action)
            {
                RHS_nSuccess = doShotNoAction(RHS_guess,RHS_solTemp,V,W0,xi,
                                   h,false_vacuum,constants,odeToSolve,events,
                                   graded_reltol_rhs,graded_abstol_rhs,thresh,
                                   overshoot_lower,overshoot_upper,bcAdjuster,
                                   fixed_background);
            }
            else
            {
                RHS_nSuccess = doShotWithAction(RHS_guess,RHS_solTemp,V,W0,xi,
                                   h,false_vacuum,constants,odeToSolveAction,
                                   eventsAction,action_reltol,action_abstol,
                                   thresh,overshoot_lower,overshoot_upper,
                                   bcAdjusterAction,
                                   fixed_background);
            }

            //Now check whether we found an overshoot or an undershoot:
            bool RHS_overshoots = check_overshoot(RHS_solTemp,RHS_upper,
                                                  RHS_lower,
                                                  RHS_guess,outStream,1,
                                                  reverse_shooting);
            std::cout << "\n" << (RHS_overshoots ? "overshoot " : "undershoot ")
                      << "\nRHS_lower = " << RHS_lower
                      << "\nRHS_upper = " << RHS_upper
                      << "\nRHS_guess = " << RHS_guess;
            if(RHS_overshoots)
            {
                overshoots_rhs_found++;
            }
            else
            {
                undershoots_rhs_found++;
            }
            printEvents(RHS_solTemp);
            end_type_rhs.push_back(RHS_solTemp.IE[RHS_solTemp.IE.size() - 1]);
            y0_rhs_list.push_back(RHS_guess);
            found_solution_rhs = (RHS_diff < tol_RHS) && compute_action
                                    && (end_type_rhs.back() == 2);
        }

        //Record locations of last undershoot and overshoot:
        if(!reverse_shooting)
        {
            last_undershoot_lhs = LHS_upper;
            last_undershoot_rhs = RHS_lower;
            last_overshoot_lhs = LHS_lower;
            last_overshoot_rhs = RHS_upper;
        }
        else
        {
            last_undershoot_lhs = LHS_lower;
            last_undershoot_rhs = RHS_upper;
            last_overshoot_lhs = LHS_upper;
            last_overshoot_rhs = RHS_lower;
        }

        //Determine matching conditions for the scale factor. This helps us
        //spot when we have converged to a pair of solutions from either side
        //that don't actually match up:
        /*
        int nEvent_ap_zero_LHS;
        int nEvent_ap_zero_RHS;
        bool LHShasap0 = false;
        bool RHShasap0 = false;
        for(int i = 0; i < LHS_solTemp.IE.size();i++)
        {
            LHShasap0 = (LHS_solTemp.IE[i] == 1);
            if(LHShasap0)
            {
                nEvent_ap_zero_LHS = i;
                break;
            }
        }
        for(int i = 0; i < RHS_solTemp.IE.size();i++)
        {
            RHShasap0 = (RHS_solTemp.IE[i] == 1);
            if(RHShasap0)
            {
                nEvent_ap_zero_RHS = i;
                break;
            }
        }
        if(LHShasap0 && RHShasap0)
        {
            value_type a_diff = LHS_solTemp.YE[nEvent_ap_zero_LHS][2]
                             - RHS_solTemp.YE[nEvent_ap_zero_RHS][2];
            std::cout << "\na_diff = " << a_diff;
            value_type a_av = (LHS_solTemp.YE[nEvent_ap_zero_LHS][2]
                               + RHS_solTemp.YE[nEvent_ap_zero_RHS][2])/_2p0;
            std::cout << "\na_diff/a_av = " << a_diff/a_av;
            a_match = abs(a_diff/(a_av)) < precision;
        }
        else
        {
            a_match = false;
        }
        */



        //Adjust bounds so that the events function can integrate more
        //efficiently, since if we pass these bounds at any point
        //with y' < 0, we
        //have definitely overshot and don't need to track a solution which
        //diverges right up until a = 0
        /*
        if(!reverse_shooting)
        {
            eventsTrack.y0upperBound = RHS_upper;
            eventsTrack.y0lowerBound = LHS_lower;
            eventsNoTrack.y0upperBound = RHS_upper;
            eventsNoTrack.y0lowerBound = LHS_lower;
        }*/
        counter++;
        /*
        outStream << "\nn = " << counter << "\nLHS "
                << (LHS_overshoots ? " overshoot " : " undershoot ")
                << "\nRHS " <<(RHS_overshoots ? " overshoot " : " undershoot ");
        */
        LHS_diff = abs(LHS_upper - LHS_lower);
        RHS_diff = abs(RHS_upper - RHS_lower);
        if(counter > nMax)
        {
            break;
        }

        //Clear solution grids for the next pass, unless we are done with
        //that half:
        if(LHS_diff > tol_LHS)
        {
            LHS_solTemp.erase();
        }
        if(RHS_diff > tol_RHS)
        {
            RHS_solTemp.erase();
        }

        //If we appear to be within tolerance, upgrade to maximum precision
        //and check consistency. If we are not consistent, then continue
        //until we are:
        if(!(LHS_diff > tol_LHS || RHS_diff > tol_RHS))
        {
            //Find consistent bounds for the LHS:
            if(inconsistent_lhs)
            {
                graded_reltol_lhs = action_reltol;
                graded_abstol_lhs = action_abstol;
                checkConsistency(y0_lhs_list,end_type_lhs,
                                 last_undershoot_lhs,last_overshoot_lhs,
                                 LHS_solTemp,V,W0,xi,h,false_vacuum,
                                 constants,odeToSolve,events,
                                 graded_reltol_lhs,graded_abstol_lhs,thresh,
                                 overshoot_lower,overshoot_upper,
                                 bcAdjuster,LHS_upper,LHS_lower,0,
                                 compute_action,fixed_background);

                if(!reverse_shooting)
                {
                    LHS_upper = last_undershoot_lhs;
                    LHS_lower = last_overshoot_lhs;
                    LHS_diff = abs(LHS_upper - LHS_lower);
                }
                else
                {
                    LHS_upper = last_overshoot_lhs;
                    LHS_lower = last_undershoot_lhs;
                    LHS_diff = abs(LHS_upper - LHS_lower);
                }
                inconsistent_lhs = false;
                LHS_diff = abs(LHS_upper - LHS_lower);
            }
            //Find consistent bounds for the RHS:
            if(inconsistent_rhs)
            {
                graded_reltol_rhs = action_reltol;
                graded_abstol_rhs = action_abstol;
                checkConsistency(y0_rhs_list,end_type_rhs,
                                 last_undershoot_rhs,last_overshoot_rhs,
                                 RHS_solTemp,V,W0,xi,h,false_vacuum,
                                 constants,odeToSolve,events,
                                 graded_reltol_rhs,graded_abstol_rhs,thresh,
                                 overshoot_lower,overshoot_upper,
                                 bcAdjuster,RHS_upper,RHS_lower,1,
                                 compute_action,fixed_background);

                if(!reverse_shooting)
                {
                    RHS_upper = last_overshoot_rhs;
                    RHS_lower = last_undershoot_rhs;
                    RHS_diff = abs(RHS_upper - RHS_lower);
                }
                else
                {
                    RHS_upper = last_undershoot_rhs;
                    RHS_lower = last_overshoot_rhs;
                    RHS_diff = abs(RHS_upper - RHS_lower);
                }
                inconsistent_rhs = false;
                RHS_diff = abs(RHS_upper - RHS_lower);
            }
            inconsistent = inconsistent_rhs || inconsistent_lhs;
            compute_action = !inconsistent;
        }
        //Make sure we terminate only on an undershoot:
        both_undershoots = (end_type_lhs.back() == 2) &&
                                    (end_type_rhs.back() == 2);
        no_undershoots_found_lhs = no_undershoots_found_lhs &&
                                    end_type_lhs.back() != 2;
        no_undershoots_found_rhs = no_undershoots_found_rhs &&
                                    end_type_rhs.back() != 2;
        found_solution = (found_solution_lhs && found_solution_rhs)
                                && !inconsistent;
        //Emergency exit:
        if(counter > 200)
        {
            break;
        }
        //Exit if we found two undershoot solution of the requires precision.
        //But don't impose this condition if we never found ANY undershoots
        //(this can occur for the Hawking-Moss solution, for example)
    }
    while(!(found_solution && (both_undershoots ||
                               (no_undershoots_found_lhs &&
                                no_undershoots_found_rhs) )));
    /*while( ((LHS_diff > tol_LHS || RHS_diff > tol_RHS) || inconsistent )
           || !both_undershoots );*/


    //We found the bounce. Are we within precision?
    if(!found_solution)
    {
        outStream << "Loop exited before finding bounce." << std::endl;
        //reportToCaller("Loop exited before finding bounce.\n");
    }




    //Bounce boundary conditions:
    /*
    dSSolTypeAction LHS_y0Action;
    bool undershoot_lhs = false;
    bool undershoot_rhs = false;
    int counter_lhs = 0;
    int counter_rhs = 0;
    std::cout.precision(50);
    initStepMax = time_type(0.01)/chimax;
    for(int i = end_type_lhs.size() - 1;i >= 0;i--)
    {
        //Search for the last undershoot we found
        if(end_type_lhs[i] != 2)
        {
            continue;
        }
        LHS_guess = y0_lhs_list[i];
        //value_type LHS_y0Action_List[] = {LHS_guess,_0p0,_0p0,_0p0,_0p0};
        value_type LHS_y0Action_List[] = {LHS_guess,_0p0,_0p0,_0p0,_0p0,_0p0};
        //LHS_y0Action.assign(LHS_y0Action_List,LHS_y0Action_List + 5);
        LHS_y0Action.assign(LHS_y0Action_List,LHS_y0Action_List + 6);
        dSSolTypeAction LHS_y0Action1step = LHS_y0Action;
        time_type LHS_epsAction;
        std::cout << "\nY_{lhs} = ";
        for(int i = 0;i < 5;i++)
        {
            std::cout << "\nY[" << i << "] = " << LHS_y0Action_List[i];
        }
        //Adjust boundary conditions to avoid co-ordinate singularities:
        //if(false)
        if(abs(V.d(LHS_guess)) < stepError && linear_BCs)
        {
            LHS_epsAction = dynBCvacAction(LHS_y0Action,LHS_y0Action1step);
        }
        else
        {
            LHS_epsAction = (*bcAdjusterAction)(LHS_y0Action,LHS_y0Action1step);
        }
        //Clear solution grids for the final result:
        LHS_solTemp.erase();
        //Time span to integrate over:
        time_type LHS_tspanArray[2] = { LHS_epsAction , chimax };
        std::vector< time_type > LHS_tspan;
        LHS_tspan.assign(LHS_tspanArray,LHS_tspanArray + 2);
        //LHS solution:
        std::cout << "\nComputing bounce solution.";
        std::cout << "\naction_reltol = " << action_reltol;
        std::cout << "\naction_abstol = " << action_abstol;
        int LHS_nSuccess;
        if(use_analytic_a)
        {
            LHS_nSuccess = ode_solve_trans<time_type,dSSolType,value_type>
                  (*odeToSolveAction,LHS_y0Action1step,LHS_tspan,*eventsAction,
                   odeSolverToUse,LHS_solTemp,action_reltol,action_abstol,
                   initStepMax,thresh,overshoot_lower,overshoot_upper,V,W0,true,
                   xi,h,false,true,false_vacuum);
        }
        else
        {
            LHS_nSuccess = ode_solve<time_type,dSSolTypeAction,value_type>
                (*odeToSolveAction,LHS_y0Action1step,LHS_tspan,*eventsAction,
                 odeSolverToUse,LHS_solTemp,action_reltol,action_abstol,
                 initStepMax);
        }
        if(LHS_nSuccess != 1)
        {
            outStream << "nSuccess = " << LHS_nSuccess
                      << " when solving for AdS-Flat"
                      << "instantons." << std::endl;
            throw "Integration failure.";
        }
        std::cout << "\nComputed " << LHS_solTemp.T.size()
                  << " points on the LHS.";
        std::cout << "\nEvents = (";
        for(int i = 0;i < LHS_solTemp.IE.size();i++)
        {
            if(i == 0)
            {
                std::cout << LHS_solTemp.IE[i];
            }
            else
            {
                std::cout << "," << LHS_solTemp.IE[i];
            }
        }
        std::cout << ").\n";
        if(LHS_solTemp.IE[LHS_solTemp.IE.size() - 1] != 2)
        {
            std::cout << "\nSolution not an undershoot at higher precision."
                    << " Keep searching.";
            continue;
        }
        else
        {
            //Found an undershoot:
            break;
        }
    }*/
    //Export LHS solution
    solOutLHS = LHS_solTemp;


    //RHS:
    /*
    for(int i = end_type_rhs.size() - 1;i >= 0;i--)
    {
        //Search for the last undershoot we found
        if(end_type_rhs[i] != 2)
        {
            continue;
        }
        RHS_guess = y0_rhs_list[i];
        //Setup BCs:
        dSSolTypeAction RHS_y0Action;
        //value_type RHS_y0Action_List[] = { RHS_guess,_0p0,_0p0,_0p0,_0p0};
        value_type RHS_y0Action_List[] = { RHS_guess,_0p0,_0p0,_0p0,_0p0,_0p0};
        //RHS_y0Action.assign(RHS_y0Action_List,RHS_y0Action_List + 5);
        RHS_y0Action.assign(RHS_y0Action_List,RHS_y0Action_List + 6);
        dSSolTypeAction RHS_y0Action1step = RHS_y0Action;
        //Output initial conditions:
        std::cout << "\nY_{rhs} = ";
        for(int i = 0;i < RHS_y0Action.size();i++)
        {
            std::cout << "\nY[" << i << "] = " << RHS_y0Action_List[i];
        }
        time_type RHS_epsAction;
        //if(false)
        if(abs(V.d(RHS_guess)) < stepError && linear_BCs)
        {
            RHS_epsAction = dynBCvacAction(RHS_y0Action,RHS_y0Action1step);
        }
        else
        {
            RHS_epsAction = (*bcAdjusterAction)(RHS_y0Action,RHS_y0Action1step);
        }
        //time_type LHS_epsAction = (*bcAdjusterAction)(LHS_y0Action,
        //                                              LHS_y0Action1step);
        //time_type RHS_epsAction = (*bcAdjusterAction)(RHS_y0Action,
        //                                              RHS_y0Action1step);

        //Clear solution grids for the final result:
        RHS_solTemp.erase();

        //Time span:
        time_type RHS_tspanArray[2] = { RHS_epsAction , chimax };
        std::vector< time_type > RHS_tspan;
        RHS_tspan.assign(RHS_tspanArray,RHS_tspanArray + 2);
        initStepMax = time_type(0.01)/chimax;

        //RHS solution
        int RHS_nSuccess;
        if(use_analytic_a)
        {
            RHS_nSuccess = ode_solve_trans<time_type,dSSolType,value_type>
                  (*odeToSolveAction,RHS_y0Action1step,RHS_tspan,*eventsAction,
                   odeSolverToUse,RHS_solTemp,action_reltol,action_abstol,
                   initStepMax,thresh,overshoot_lower,overshoot_upper,V,W0,true,
                   xi,h,true,true,false_vacuum);
        }
        else
        {
            RHS_nSuccess = ode_solve<time_type,dSSolTypeAction,value_type>
                (*odeToSolveAction,RHS_y0Action1step,RHS_tspan,*eventsAction,
                 odeSolverToUse,RHS_solTemp,action_reltol,action_abstol,
                 initStepMax);
        }
        if(RHS_nSuccess != 1)
        {
            outStream << "nSuccess = " << RHS_nSuccess
                      << " when solving for AdS-Flat"
                      << "instantons." << std::endl;
            throw "Integration failure.";
        }
        std::cout << "\nComputed " << RHS_solTemp.T.size()
                  << " points on the RHS.";
        std::cout << "\nEvents = (";
        for(int i = 0;i < RHS_solTemp.IE.size();i++)
        {
            if(i == 0)
            {
                std::cout << RHS_solTemp.IE[i];
            }
            else
            {
                std::cout << "," << RHS_solTemp.IE[i];
            }
        }
        std::cout << ").\n";
        //Check if this is a undershoot:
        if(RHS_solTemp.IE[RHS_solTemp.IE.size() - 1] != 2)
        {
            std::cout << "\nSolution not an undershoot at higher precision."
                    << " Keep searching.";
            continue;
        }
        else
        {
            //Found an undershoot:
            break;
        }
    }*/
    //Export RHS solution
    solOutRHS = RHS_solTemp;

    //LHS_guess = (LHS_upper + LHS_lower)/_2p0;
    //RHS_guess = (RHS_upper + RHS_lower)/_2p0;
    //value_type LHS_y0Action_List[] = { LHS_guess , _0p0 , _0p0 ,a_scale,_0p0};


    //value_type RHS_y0Action_List[] = { RHS_guess , _0p0 , _0p0 ,a_scale,_0p0};











    //solution_grid< time_type , dSSolTypeAction , value_type > LHS_solTemp;
    //solution_grid< time_type , dSSolTypeAction , value_type > RHS_solTemp;

    //std::cout << "\nSetup bounce time-span.";

    //Time span to integrate over:















    //Splice solutions.
    //First copy in the RHS solution up to the point where a' = 0.

    //Retrieve the point where a' = 0:
    int nEvent_ap_zero;
    bool hasap0 = false;
    for(int i = 0; i < RHS_solTemp.IE.size();i++)
    {
        hasap0 = (RHS_solTemp.IE[i] == 1);
        if(hasap0)
        {
            nEvent_ap_zero = i;
            break;
        }
    }
    if(!hasap0)
    {
        outStream << "\nFailed to find matching point for RHS solution."
                  << "Export RHS solution.";
        solOut.Y = RHS_solTemp.Y;
        solOut.T = RHS_solTemp.T;
        solOut.YE = RHS_solTemp.YE;
        solOut.TE = RHS_solTemp.TE;
        solOut.IE = RHS_solTemp.IE;
        DSout = !fixed_background ? solOut.Y[solOut.Y.size() - 1][4] :
             solOut.Y[solOut.Y.size() - 1][2];
        throw "RHS solution does not reach a' = 0.";
    }

    //False vacuum action.
    value_type yfv2 = false_vacuum*false_vacuum;
    value_type h2 = h*h;
    value_type factor = (_1p0 - xi*h2*yfv2);
    value_type factor2 = _1p0 - xi*(_1p0 + _6p0*xi)*h2*yfv2;
    value_type factor3 = _1p0 - xi*(_1p0 - _6p0*xi)*h2*yfv2;
    //value_type S0 = -_24p0*PI*PI*factor*factor*factor2/(factor3*h2*h2*W0);
    //For the moment, this only works with xi = 0:
    //Determine the total time period:

    //Retrieve the point where a' = 0:
    int LHS_nEvent_ap_zero;
    hasap0 = false;
    for(int i = 0; i < LHS_solTemp.IE.size();i++)
    {
        hasap0 = (LHS_solTemp.IE[i] == 1);
        if(hasap0)
        {
            LHS_nEvent_ap_zero = i;
            break;
        }
    }
    std::cout << "\nOk here! 2.";
    if(!hasap0)
    {
        outStream << "\nFailed to find matching point for LHS solution."
                  << "Export LHS solution.";
        solOut.Y = LHS_solTemp.Y;
        solOut.T = LHS_solTemp.T;
        solOut.YE = LHS_solTemp.YE;
        solOut.TE = LHS_solTemp.TE;
        solOut.IE = LHS_solTemp.IE;
        DSout = !fixed_background ? solOut.Y[solOut.Y.size() - 1][4] :
                solOut.Y[solOut.Y.size() - 1][2];
        throw "LHS solution does not reach a' = 0.";
    }
    std::cout << "\nOk here! 3.";
    time_type chimid = RHS_solTemp.TE[nEvent_ap_zero];
    time_type chimax_m_chimid = LHS_solTemp.TE[LHS_nEvent_ap_zero];
    time_type chi_m = !fixed_background ? chimid + chimax_m_chimid : PI/H0;

    //value_type H0 = h*sqrt(W0/_3p0);
    /*
    value_type chi_max = value_type(chi_m);
    value_type chimid = value_type(RHS_solTemp.TE[nEvent_ap_zero]);
    value_type X1 = cos(H0*(chi_max - chimid));
    value_type X2 = cos(H0*chimid);
    value_type root3 = sqrt(_3p0);
    value_type S0 = _6p0*PI*PI*( X1*(X1 + root3)*(X1 - root3) +
                                  X2*(X2 + root3)*(X2 - root3) )/W0;
    */
    value_type chi_max = value_type(chi_m);
    value_type B1_rhs,B2_rhs,B1_lhs,B2_lhs,B1_rhs_half,B1_lhs_half,B3;
    value_type cos_factor = cos(H0*chi_max);
    value_type cos_ratio = (_1p0 + cos_factor)/H0;
    if(!fixed_background)
    {
        B1_rhs = RHS_solTemp.YE[RHS_solTemp.YE.size() - 1][4];
        B2_rhs = RHS_solTemp.YE[RHS_solTemp.YE.size() - 1][5];
        B1_lhs = LHS_solTemp.YE[LHS_solTemp.YE.size() - 1][4];
        B2_lhs = LHS_solTemp.YE[LHS_solTemp.YE.size() - 1][5];
        B1_rhs_half = RHS_solTemp.YE[nEvent_ap_zero][4];
        B1_lhs_half = LHS_solTemp.YE[LHS_nEvent_ap_zero][4];
        B3 = -_2p0*PI*PI*h*h*cos_ratio*cos_ratio*(cos_factor - _2p0);
    }
    else
    {
        B1_rhs = RHS_solTemp.YE[RHS_solTemp.YE.size() - 1][2];
        B2_rhs = _0p0;
        B1_lhs = LHS_solTemp.YE[LHS_solTemp.YE.size() - 1][2];
        B2_lhs = _0p0;
        B1_rhs_half = RHS_solTemp.YE[nEvent_ap_zero][2];
        B1_lhs_half = LHS_solTemp.YE[LHS_nEvent_ap_zero][2];
        B3 = _0p0;
    }
    //std::cout << "\nS0 = " << S0;


    //Get everything before the point where a' = 0:
    solOut.terminated_early = RHS_solTemp.terminated_early;
    //solOut.Y.reserve(RHS_solTemp.eventLocations[nEvent_ap_zero]);
    //solOut.T.reserve(RHS_solTemp.eventLocations[nEvent_ap_zero]);
    int nLHS = 0;
    value_type a_scale_3 = a_scale*a_scale*a_scale;
    if(!fixed_background)
    {
        for(int i = 0;i < RHS_solTemp.Y.size();i++)
        {
            /*
            if(int(RHS_solTemp.Y[i].size()) != 5)
            {
                std::cout << "\nError - size not equal to 5!";
            }
            std::cout.precision(5);
            for(int j = 0;j < RHS_solTemp.Y[i].size();j++)
            {
                std::cout << "\ny[" << j << "] = " << RHS_solTemp.Y[i][j];
            }*/
            //Get everything before we encounter a' = 0:
            if(RHS_solTemp.Y[i][3] + cos(H0*value_type(RHS_solTemp.T[i]))> _0p0)
            {
                std::vector<value_type> vec = RHS_solTemp.Y[i];
                value_type t = value_type(RHS_solTemp.T[i]);
                //Modify delta_a so that it matches in the middle (hopefully):
                value_type H0chimid = H0*value_type(chimid);
                value_type H0chimax = H0*chi_max;
                RHS_solTemp.Y[i][2] += (sin(H0*t)
                            - sin(H0*(chi_max - t)))/H0;
                RHS_solTemp.Y[i][3] += (cos(H0*t)
                            + cos(H0*(chi_max - t)));
                RHS_solTemp.Y[i][2] /= a_scale;
                RHS_solTemp.Y[i][3] /= a_scale;
                //std::cout.precision(50);
                RHS_solTemp.Y[i][4] /= a_scale_3;
                RHS_solTemp.Y[i][5] /= a_scale_3;
                //Modify component 5 so that it will match the LHS version:
                RHS_solTemp.Y[i][5] += (sin(H0*t)
                            - sin(H0*(chi_max - t)) )/H0;
                solOut.Y.push_back(RHS_solTemp.Y[i]);

                //solOut.Y[i][4] /= a_scale_3;//Counter-act the
                //effect of changing
                    //a in the action.
                    /*
                if(!use_analytic_a)
                {
                    solOut.Y[i][4] += - S0;//Subtract off false
                                        //vacuum action and include the other
                                        //correction.
                }*/
                solOut.T.push_back(t);
                nLHS ++;
            }
            else
            {
                break;
            }
        }
    }
    else
    {
        for(int i = 0;i < RHS_solTemp.Y.size();i++)
        {
            //Get everything before we encounter a' = 0:
            if(cos(H0*value_type(RHS_solTemp.T[i]))> _0p0)
            {
                std::vector<value_type> vec = RHS_solTemp.Y[i];
                value_type t = value_type(RHS_solTemp.T[i]);
                //Modify delta_a so that it matches in the middle (hopefully):
                value_type H0chimid = H0*value_type(chimid);
                value_type H0chimax = H0*chi_max;
                RHS_solTemp.Y[i][2] /= a_scale_3;
                solOut.Y.push_back(RHS_solTemp.Y[i]);
                solOut.T.push_back(t);
                nLHS ++;
            }
            else
            {
                break;
            }
        }
    }
    //Copy all the events, because a' = 0 is terminal anyway:
    solOut.YE = RHS_solTemp.YE;
    solOut.TE = RHS_solTemp.TE;
    solOut.IE = RHS_solTemp.IE;

    std::cout << "\nOk here! 1.";

    //Test integration of B2 using the complete solution, via trapezium
    //rule estimate:
    //STILL WORK TO BE DONE MODIFYING THIS FOR COMPATIBILITY WITH
    //FIXED BACKGROUND CASE
    std::vector<value_type> integrand (int(solOut.T.size()),_0p0);
    std::vector<value_type> chi_list (int(solOut.T.size()),_0p0);
    std::vector<value_type> action_integrand (int(solOut.T.size()),_0p0);
    value_type prefactor = -_6p0*PI*PI*h*h;
    value_type prefactor2 = -_2p0*PI*PI;
    value_type B2_estimate,B1_estimate;
    if(!fixed_background)
    {
        for(int i = 0;i < integrand.size();i++)
        {
            value_type sinarg = sin(H0*value_type(solOut.T[i]));
            value_type delta_a = solOut.Y[i][2];
            value_type a = delta_a + sin(H0*(chi_max -
                                             value_type(solOut.T[i])))/H0;

            integrand[i] = prefactor*(_3p0*sinarg*sinarg*delta_a +
                                      _3p0*H0*sinarg*delta_a*delta_a +
                                      H0*H0*delta_a*delta_a*delta_a);
            chi_list[i] = value_type(solOut.T[i]);
            action_integrand[i] = prefactor2*a*a*a*V(solOut.Y[i][0]);
        }
        B2_estimate = trap(chi_list,integrand);
        B1_estimate = trap(chi_list,action_integrand);
    }
    else
    {
        B2_estimate = _0p0;
        for(int i = 0;i < integrand.size();i++)
        {
            value_type sinarg = sin(H0*value_type(solOut.T[i]));
            value_type a = sin(H0*(chi_max -
                                             value_type(solOut.T[i])))/H0;
            chi_list[i] = value_type(solOut.T[i]);
            action_integrand[i] = prefactor2*a*a*a*V(solOut.Y[i][0]);
        }
        B1_estimate = trap(chi_list,action_integrand);
    }
    std::cout << "\nTrapezium rule estimate of B2: B2(H0) = "
              << B2_estimate;
    std::cout << "\nTrapezium rule estimate of B1: B1(H0) = "
              << B1_estimate;

    //Now comes the trickier part - we add on the LHS solution. However, we have
    // to reverse the order, flipping the signs of the derivatives, and shifting
    //the time variable up, ie, chi -> chi_m - chi where chi_m is the total
    //length of the combined solution, in other words, the time for the RHS
    //solution to reach a' = 0 plus the time for the LHS solution to reach
    //a' = 0.



    //The difference in the actions at the centre represents the true action of
    //the solution, since we start both sides with S = 0, but we can always
    //add an arbitrary constant to force the solutions to match in the middle.
    //Note that by definition in the equations we solve, this takes into account
    //the difference between the action of the bounce and that of the false
    //vacuum.
    /*
    DSout = RHS_solTemp.YE[nEvent_ap_zero][4] +
            LHS_solTemp.YE[LHS_nEvent_ap_zero][4] - S0;*/
    //Best guess for the decay exponent uses half of B1 from the RHS, half from
    //the lhs, and the whole of B2 from the LHS (the RHS is unreliable, as
    //delta_a does not vanish at chi_max there). B3 is determined analytically
    //and B = B1 + B2 + B3.
    std::cout << "\nB1 = " << B1_rhs_half + B1_lhs_half;
    std::cout << "\nB2 = " << B2_lhs;
    std::cout << "\nB3 = " << B3;
    DSout = B1_rhs_half + B1_lhs_half + B2_lhs + B3;
    if(!fixed_background)
    {
        std::cout << "\nRHS_contrib = " << RHS_solTemp.YE[nEvent_ap_zero][4];
        std::cout << "\nLHS_contrib = "
                  << LHS_solTemp.YE[LHS_nEvent_ap_zero][4];
    }
    else
    {
        std::cout << "\nRHS_contrib = " << RHS_solTemp.YE[nEvent_ap_zero][2];
        std::cout << "\nLHS_contrib = "
                  << LHS_solTemp.YE[LHS_nEvent_ap_zero][2];
    }

    if(use_analytic_a)
    {
        DSout += RHS_solTemp.S_correction
                                + LHS_solTemp.S_correction;
        std::cout << "\nRHS_solTemp.S_correction = "
                  << RHS_solTemp.S_correction;
        std::cout << "\nLHS_solTemp.S_correction = "
                  << LHS_solTemp.S_correction;
    }
    std::cout << "\nOk here! 4.";



    //Add the RHS solutions, flipping as necessary, up to but NOT including a'=0
    /*solOut.Y.reserve(RHS_solTemp.eventLocations[nEvent_ap_zero] +
                     LHS_solTemp.eventLocations[LHS_nEvent_ap_zero]);
    solOut.T.reserve(RHS_solTemp.eventLocations[nEvent_ap_zero] +
                     LHS_solTemp.eventLocations[LHS_nEvent_ap_zero]);*/
    //First find the location of the last point before a' = 0.
    int nRHS = 0;
    bool rhs_condition;
    for(int i = 0;i < LHS_solTemp.Y.size();i++)
    {
        rhs_condition = !fixed_background ? (LHS_solTemp.Y[i][3] +
                        cos(H0*value_type(LHS_solTemp.T[i])) > _0p0) :
                        cos(H0*value_type(LHS_solTemp.T[i])) > _0p0;
        if(rhs_condition)
        {
            nRHS++;
        }
        else
        {
            break;
        }
    }
    std::cout << "\nOk here! 5.";
    for(int i = nRHS - 1;i >= 0;i--)
    {
        /*
        if(int(LHS_solTemp.Y[i].size()) != 5)
        {
            std::cout << "\nError - size not equal to 5!";
            std::cout << "\ni = " << i;
            std::cout << "\nSize = " << int(LHS_solTemp.Y[i].size());
        }
        std::cout.precision(5);
        for(int j = 0;j < LHS_solTemp.Y[i].size();j++)
        {
            std::cout << "\ny[" << j << "] = " << RHS_solTemp.Y[i][j];
        }*/
        //Flip derivative signs as we are going in the opposite direction:
        if(!fixed_background)
        {
            LHS_solTemp.Y[i][1] *= -_1p0;
            LHS_solTemp.Y[i][3] *= -_1p0;
            LHS_solTemp.Y[i][4] *= -_1p0;
            LHS_solTemp.Y[i][5] *= -_1p0;

            //Modify a:
            LHS_solTemp.Y[i][2] /= a_scale;
            LHS_solTemp.Y[i][3] /= a_scale;
            LHS_solTemp.Y[i][4] /= a_scale_3;
            LHS_solTemp.Y[i][5] /= a_scale_3;
        }
        else
        {
            LHS_solTemp.Y[i][1] *= -_1p0;

            //Modify a:
            LHS_solTemp.Y[i][2] /= a_scale_3;
        }

        //Shift the action up to match the RHS solution in the middle:
        /*
        if(use_analytic_a)
        {
            LHS_solTemp.Y[i][4] += DSout + RHS_solTemp.S_correction
                                + LHS_solTemp.S_correction;
        }
        else
        {

        }
        */
        if(!use_analytic_a)
        {
            //LHS_solTemp.Y[i][4] += DSout;
            if(!fixed_background)
            {
                LHS_solTemp.Y[i][4] += B1_rhs_half + B1_lhs_half;
            }
            else
            {
                LHS_solTemp.Y[i][2] += B1_rhs_half + B1_lhs_half;
            }
        }
            //Note that this includes the subtraction of the false
            //vacuum action.




        //Copy solution across, with the adjusted time variable:
        if(!fixed_background)
        {
            solOut.Y.push_back(LHS_solTemp.Y[i]);
            solOut.T.push_back(chi_m - LHS_solTemp.T[i]);
        }
        else
        {
            time_type t_i = chi_m - LHS_solTemp.T[i];
            solOut.T.push_back(t_i);
            std::vector<value_type> value_i (5,_0p0);
            value_i[0] = LHS_solTemp.Y[i][0];
            value_i[1] = LHS_solTemp.Y[i][1];
            value_i[2] = sin(H0*value_type(t_i))/H0;
            value_i[3] = cos(H0*value_type(t_i));
            value_i[4] = LHS_solTemp.Y[i][2];
            solOut.Y.push_back(value_i);
        }
    }
    std::cout << "\nOk here! 6.";

    //Number of elements selected from each solution:
    //int nLHS = LHS_solTemp.eventLocations[LHS_nEvent_ap_zero];
    //int nRHS = RHS_solTemp.eventLocations[nEvent_ap_zero] + 1;//+1 as a' =0
        //included in this solution.
    solOut.nSteps = nLHS + nRHS - 1;
    //Copy events functions in the same way as the solutions:
    /*solOut.YE.reserve(nEvent_ap_zero + LHS_nEvent_ap_zero);
    solOut.TE.reserve(nEvent_ap_zero + LHS_nEvent_ap_zero);
    solOut.IE.reserve(nEvent_ap_zero + LHS_nEvent_ap_zero);*/
    //solOut.eventLocations.reserve(nEvent_ap_zero + LHS_nEvent_ap_zero);
    for(int i = int(LHS_solTemp.YE.size()) - 1; i >= 0;i--)
    {
        std::cout << "\nLHS_solTemp.YE[" << i << "].size() = "
                  << int(LHS_solTemp.YE[i].size());
        LHS_solTemp.YE[i][1] *= -_1p0;
        if(!fixed_background)
        {
            LHS_solTemp.YE[i][3] *= -_1p0;
        }
        /*
        if(!use_analytic_a)
        {
            LHS_solTemp.YE[i][4] += RHS_solTemp.YE[nEvent_ap_zero][4] +
                        LHS_solTemp.YE[LHS_nEvent_ap_zero][4] - S0;
        }*/
        time_type t_i = chi_m - LHS_solTemp.TE[i];
        std::vector<value_type> value_i (5,_0p0);
        value_i[0] = LHS_solTemp.YE[i][0];
        value_i[1] = LHS_solTemp.YE[i][1];
        value_i[2] = sin(H0*value_type(t_i))/H0;
        value_i[3] = cos(H0*value_type(t_i));
        value_i[4] = LHS_solTemp.YE[i][2];
        solOut.YE.push_back(value_i);
        solOut.TE.push_back(t_i);
        solOut.IE.push_back(LHS_solTemp.IE[i]);
        //Adjust locations of LHS events:
        /*solOut.eventLocations.push_back(nLHS + nRHS
                                        - 1 -  RHS_solTemp.eventLocations[i]);*/
    }
    std::cout << "\nOk here! 7.";

    //Return any additional return values we might want:
    //B1:
    other_return_values[0] = B1_rhs_half + B1_lhs_half;
    //B2:
    other_return_values[1] = B2_lhs;
    //B3:
    other_return_values[2] = B3;
    //B1rhs:
    other_return_values[3] = B1_rhs_half;
    //B1lhs:
    other_return_values[4] = B1_lhs_half;
    //B2 (RHS version, probably inaccurate):
    other_return_values[5] = B2_rhs;
    //B1_estimate
    other_return_values[6] = B1_estimate;
    //B2_estimate
    other_return_values[7] = B2_estimate;

    outStream << "H0*chi_max/pi - 1: "
              << H0*solOut.T[int(solOut.T.size()) - 1]/PI - _1p0;

    outStream << "Done." << std::endl;
    //Done.
    }
    catch(char* c)
    {
        outStream << c;
    }
    catch(...)
    {
        outStream << "An unknown error occurred.";
    }
}
//------------------------------------------------------------------------------
template< class value_type , class time_type >
void odeSolve_dSsingle(value_type y0, value_type false_vacuum ,
                           value_type true_vacuum , value_type barrier ,
                           potential< value_type >& V , time_type chimax ,
                           int odeSolverToUse,
                           value_type RelTol, value_type AbsTol,
                           value_type stepError, value_type xi,
                           solution_grid< time_type , std::vector< value_type >,
                           value_type >& solOut,value_type& DSout,
                           value_type precision,value_type h,value_type W0,
                           std::ostream& outStream,
                           bool track_scale_factor,bool linear_BCs,
                           bool fixed_background)
{
//INPUTS:



    try
    {
     outStream << "Starting Integration for y0 = \n" << y0 << std::endl;
     outStream << "\nV(y0) = " << V(y0) << std::endl;

    //Ok. we proceed by bisecting until we have a solution which overshoots
    //beyond
    //chimax, which stands in for infinity.

    //Overshoots occur when the solution crosses zero before chimax.
    //Undershoots when
    //its derivative crosses zero before then.

    //Instantiate the events function, state-type, dynamic BCs etc...
    typedef std::vector< value_type > dSSolType;//solution without action,
        //for speed.
    typedef std::vector< value_type > dSSolTypeAction;

    typedef ode_dS_delta_a< dSSolType , time_type , value_type > odeRHS_dS;
    typedef events_delta_a< dSSolType , time_type , value_type > ev_dS;
    typedef dynamicBCsAdSFlat_taylor_delta_a< value_type , dSSolType,time_type >
            dynBC_dS;
    typedef ode_dS_action_delta_a< dSSolTypeAction , time_type , value_type >
            odeRHS_dSAction;
    typedef events_delta_a< dSSolTypeAction , time_type , value_type >
            ev_dSAction;
    typedef dynamicBCsAdSFlat_taylor_delta_a< value_type , dSSolTypeAction ,
                                      time_type > dynBCAdSFlatAction;

    /*
    typedef ode_dS_delta_a< dSSolType , time_type , value_type > odeRHSAdSFlat;
    typedef eventsMethod0< dSSolType , time_type , value_type > ev_dS;
    typedef dynamicBCsAdSFlat_taylor< value_type , dSSolType , time_type >
            dynBC_dS;
    typedef ode_dS_action_delta_a< dSSolTypeAction , time_type , value_type >
            odeRHS_dSAction;
    typedef eventsMethod0< dSSolTypeAction , time_type , value_type >
            ev_dSAction;
            */
    typedef odeGravFriction< dSSolType , time_type , value_type > odeGrav;
    typedef dynamicBCsGravFriction_taylor< value_type , dSSolType ,
                                      time_type > dynBCGrav;
    typedef odeGravFrictionAction<dSSolType,time_type,value_type>
                                    odeGravAction;
    typedef eventsGravFriction<dSSolType,time_type,value_type> eventsGrav;


    //Fixed background approximation ode functions:
    typedef ode_dS_fixed< dSSolType , time_type , value_type > odeRHS_dS_fixed;
    typedef ode_dS_fixed_action< dSSolType , time_type , value_type >
            odeRHS_dS_fixed_action;
    typedef events_dS_fixed< dSSolTypeAction , time_type , value_type >
            ev_dS_fixed;
    typedef dynamicBCs_dS_fixed< value_type , dSSolTypeAction ,
                                      time_type > dynBC_dS_fixed;





    //Now compute the action and splice the two solutions together to
    //form a coherent solution.
    const value_type _2p0 = value_type(2.0);
    odeRHS_dSAction odeToSolveActionTrack(V,xi,h,W0);
    const value_type _3p0 = value_type(3.0);
    value_type H0 = h*sqrt(W0/_3p0);
    ev_dSAction eventsActionTrack(false_vacuum,true_vacuum,H0);
    dynBC_dS bcAdjusterActionTrack(xi,h,W0,V,stepError,2);
    odeGravAction odeToSolveActionNoTrack(V,xi,h,W0);
    eventsGrav eventsActionNoTrack(false_vacuum,true_vacuum);
    dynBCGrav bcAdjusterActionNoTrack(xi,h,W0,V,stepError,2);

    //Fixed background approximation:
    odeRHS_dS_fixed_action odeToSolveFixedAction(V,xi, h, W0);
    ev_dS_fixed eventsActionFixed(false_vacuum,_2p0*true_vacuum,H0);
    eventsActionFixed.dynamicSwitch[0] = false;
    dynBC_dS_fixed bcAdjusterActionFixed(xi,h,W0,V,stepError,2);


    //Switch off termination for the a' = 0 event because we usually want
    //to examine the entire solution when we are doing a single shot.
    eventsActionTrack.dynamicSwitch[0] = false;
    eventsActionNoTrack.dynamicSwitch[0] = false;

    //Switch off event termination for the y' = 0 events, so that only
    //bound crossing actually terminates us:
    //eventsActionTrack.dynamicSwitch[1] = false;
    //eventsActionNoTrack.dynamicSwitch[1] = false;


    dynamicBCs< value_type , std::vector<value_type> , time_type >*
        bcAdjusterAction;
    dynamicBCs_basicMethod< value_type, dSSolType, time_type >
        dynBCvac(xi,h,W0,V,stepError);
    dynamicBCs_actionMethod< value_type, dSSolTypeAction, time_type >
        dynBCvacAction(xi,h,W0,V,stepError);
    events_function< std::vector<value_type> , time_type, value_type >*
        eventsAction;
    odeRHS<std::vector<value_type>,time_type,value_type>* odeToSolveAction;

    if(fixed_background)
    {
        odeToSolveAction = &odeToSolveFixedAction;
        bcAdjusterAction = &bcAdjusterActionFixed;
        eventsAction = &eventsActionFixed;
    }
    else if(track_scale_factor)
    {
        odeToSolveAction = &odeToSolveActionTrack;
        bcAdjusterAction = &bcAdjusterActionTrack;
        eventsAction = &eventsActionTrack;
    }
    else
    {
        odeToSolveAction = &odeToSolveActionNoTrack;
        bcAdjusterAction = &bcAdjusterActionNoTrack;
        eventsAction = &eventsActionNoTrack;
    }

    //Setup bounds:



    //Instantiation for overshoot/undershoot:

    //Compute solution:

    //Setup ode:
    //odeRHS_dSAction odeToSolveAction(V,xi,h,W0);
    //Setup function to detect events:
    //ev_dSAction (*eventsAction)(false_vacuum,true_vacuum);

    //Setup function which translates the initial conditions at chi = 0
    //(where numerics will fail) into initial conditions at a point slightly
    //displaced from chi = 0, using analytic means:
    //dynBCAdSFlatAction bcAdjusterAction(xi,h,W0,V,stepError,2);

    bool use_analytic_a = false;
    time_type initStepMax = time_type(0.01)/chimax;
    bool reverse_shooting = false;
    shooting_constants<value_type,time_type> constants
    (stepError,V,chimax,linear_BCs, use_analytic_a,odeSolverToUse,initStepMax,
     reverse_shooting);

    int nSuccess = doShotWithAction(y0,solOut,V,W0,xi,
                                   h,false_vacuum,constants,odeToSolveAction,
                                   eventsAction,RelTol,AbsTol,
                                   RelTol,false_vacuum,true_vacuum,
                                   bcAdjusterAction,fixed_background);

    //Define a few useful constants:
    const value_type _0p0 = value_type(0.0);
    const value_type _1p0 = value_type(1.0);
    const time_type _0p01 = time_type(0.01);


    /*
    //Initial conditions for integration:
    dSSolTypeAction y0Action;
    //value_type y0ActionList[] = { y0 , _0p0 , _0p0 , _1p0 , _0p0 };
    //value_type y0ActionList[] = { y0 , _0p0 , _0p0 , _0p0 , _0p0 };
    value_type y0ActionList[] = { y0 , _0p0 , _0p0 , _0p0 , _0p0 , _0p0};
    y0Action.assign(y0ActionList,y0ActionList + 6);
    outStream << "\ny0 = " << y0;
    outStream << "\nfv = " << false_vacuum;
    //Adjust the initial conditions using bcAdjusterAction:
    dSSolTypeAction y0Action1step = y0Action;
    outStream << "Compute boundary conditions.\n";
    time_type epsAction;
    //if(abs(V.d(y0)) < stepError)
    //if(false)
    if(abs(V.d(y0)) < stepError && linear_BCs)
    {
        epsAction = dynBCvacAction(y0Action,y0Action1step);
        for(int i = 0;i < y0Action1step.size();i++)
        {
            std::cout << "\ny0Action1step[" << i << "] = " << y0Action1step[i];
        }
        for(int i = 0;i < y0Action.size();i++)
        {
            std::cout << "\ny0Action[" << i << "] = " << y0Action[i];
        }
    }
    else
    {
        epsAction = (*bcAdjusterAction)(y0Action,y0Action1step);
    }
    //time_type epsAction = (*bcAdjusterAction)(y0Action,y0Action1step);
    outStream << "Done.";
    outStream << "\ny0 = " << y0Action1step[0];
    outStream << "\ny0p = " << y0Action1step[1];
    outStream << "\na0 = " << y0Action1step[2];
    outStream << "\na0p = " << y0Action1step[3];
    outStream << "\nB1 = " << y0Action1step[4];
    outStream << "\nB2 = " << y0Action1step[5];

    //Setup range over which to integrate (assuming no termination, which will
    //almost always happen in practice)
    time_type tspanArray[2] = { epsAction , chimax };
    std::vector< time_type > tspan;
    tspan.assign(tspanArray,tspanArray + 2);
    //Maximum size of initial steps that we allow:
    time_type initStepMax = _0p01/time_type(chimax);

    //Perform the integration, saving the result in solOut.
    outStream << "Starting integration.";
	std::cout << "\nRelTol = " << RelTol;
	std::cout << "\nAbsTol = " << AbsTol;
	std::cout << "\ninitStepMax = " << initStepMax;
	std::cout << "\nchimax = " << chimax;
    int nSuccess = ode_solve<time_type,dSSolTypeAction,value_type>
    (*odeToSolveAction,y0Action1step,tspan,*eventsAction,odeSolverToUse,solOut,
     RelTol,AbsTol,initStepMax);
     */
    if(nSuccess != 1)
    {
        outStream << "nSuccess = " << nSuccess << " when solving for AdS-Flat"
                  << " instantons." << std::endl;
        throw "Integration failure.";
    }

    //Modify the action to include the contribution from the false vacuum
    //solution:
    const value_type _6p0 = value_type(6.0);
    //const value_type _3p0 = value_type(3.0);
    const value_type _24p0 = value_type(24.0);
    const value_type PI = boost::math::constants::pi<value_type>();
    value_type yfv2 = false_vacuum*false_vacuum;
    value_type h2 = h*h;
    value_type factor = (_1p0 - xi*h2*yfv2);
    value_type factor2 = _1p0 - xi*(_1p0 + _6p0*xi)*h2*yfv2;
    value_type factor3 = _1p0 - xi*(_1p0 - _6p0*xi)*h2*yfv2;
    value_type S0 = -_24p0*PI*PI*factor*factor*factor2/(factor3*h2*h2*W0);
    for(int i = 0;i < int(solOut.T.size());i++)
    {
        //solOut.Y[i][4] -= S0;
    }
    //Output the action:
    DSout = !fixed_background ? solOut.Y[int(solOut.T.size()) - 1][4] :
                                solOut.Y[int(solOut.T.size()) - 1][2];
    //Done.
    }
    catch(char* c)
    {
        outStream << c;
    }
    catch(...)
    {
        outStream << "An unknown error occurred.";
    }
}
//------------------------------------------------------------------------------
//Function which inspects a solution to determine whether it is an overshoot
//(true) or undershoot (false)
template< class value_type , class time_type >
bool check_overshoot
    (solution_grid< time_type , std::vector< value_type > ,
     value_type >& solTemp,
     value_type& upper,value_type& lower,const value_type& guess,
     std::ostream& outStream,int side,bool reverse_shoot,bool override_order)
{
    bool overshoots;
    for(int j = 0; j < int(solTemp.IE.size());j++)
    {
        outStream << "\nIE[" << j << "] = " << solTemp.IE[j];
        switch(solTemp.IE[j])
        {
        case 1:
            //a' = 0. Nothing to worry about. This is chosen as the matching
            //point. Do nothing.
            break;
        case 2:
            //undershoot, y' = 0 (note that the solution would have terminated
            // if it encountered a = 0, so finding y' = 0 indicates we encounter
            //y' = 0 before a = 0, which is the signature of an undershoot
            //solution.
            overshoots = false;
            if(side == 0)//LHS
            {
                if(!reverse_shoot)
                {
                    if(guess < upper || override_order)
                    {
                        upper = guess;//Undershoot, so need to start further
                        //from barrier.
                    }
                }
                else
                {
                    if(guess > lower || override_order)
                    {
                        lower = guess;
                    }
                }
            }
            else
            {
                if(!reverse_shoot)
                {
                    if(guess > lower || override_order)
                    {
                        lower = guess;//Undershoot, so need to start further
                        //from barrier.
                    }
                }
                else
                {
                    if(guess < upper || override_order)
                    {
                        upper = guess;
                    }
                }
            }
            break;
        case 3:
            //Found a = 0 point. Since we terminate if we find y' = 0, this
            //means we encountered a = 0 first, which is the signature of an
            //overshoot.
            overshoots = true;
            if(side == 0)
            {
                if(!reverse_shoot)
                {
                    if(guess > lower || override_order)
                    {
                        lower = guess;
                    }
                }
                else
                {
                    if(guess < upper || override_order)
                    {
                        upper = guess;
                    }
                }
            }
            else
            {
                if(!reverse_shoot)
                {
                    if(guess > upper || override_order)
                    {
                        upper = guess;
                    }
                }
                else
                {
                    if(guess < lower || override_order)
                    {
                        lower = guess;
                    }
                }
            }
        case 4:
            //Passed through the upper bound with y'<0, indicating an overshoot.
            //for the LHS solution. Thus we choose guess as the lower bound.
            overshoots = true;
            if(side == 0)
            {
                if(!reverse_shoot)
                {
                    if(guess > lower || override_order)
                    {
                        lower = guess;
                    }
                }
                else
                {
                    if(guess < upper || override_order)
                    {
                        upper = guess;
                    }
                }
            }
            else
            {
                if(!reverse_shoot)
                {
                    if(guess < upper || override_order)
                    {
                        upper = guess;
                    }
                }
                else
                {
                    if(guess > lower || override_order)
                    {
                        lower = guess;
                    }
                }
            }
            break;
        case 5:
            //Crossed lower bound with y'<0. For LHS side solutions, this
            //cannot occur without passing through y'=0 first. For RHS
            //solutions,
            overshoots = true;
            if(side == 0)
            {
                if(!reverse_shoot)
                {
                    if(guess > lower || override_order)
                    {
                        lower = guess;
                    }
                }
                else
                {
                    if(guess < upper || override_order)
                    {
                        upper = guess;
                    }
                }
            }
            else
            {
                if(!reverse_shoot)
                {
                    if(guess < upper || override_order)
                    {
                        upper = guess;
                    }
                }
                else
                {
                    if(guess > lower || override_order)
                    {
                        lower = guess;
                    }
                }
            }
            break;
        case 6:
            //Indicates a transfer to a different integration mode. Neither
            //undershoot nor overshoot.
            break;
        default:
            outStream << "Unrecognised ode-event." << std::endl;
            outStream << "Neither undershoot nor overshoot detected. "
                      << "(Perhaps integration range should be extended?) "
                      << "Relevant data:\n";
            outStream << "no. of events = " << int(solTemp.IE.size())
                      << std::endl;
            for(int i = 0; i < int(solTemp.IE.size());i++)
            {
                outStream << "IE[" << i << "] = " << solTemp.IE[i] << std::endl;
                outStream << "YE[" << i << "] = " << std::endl;
                for(int j = 0;j < solTemp.YE[i].size();j++)
                {
                    outStream << solTemp.YE[i][j] << std::endl;
                }
                outStream << "TE[" << i << "] = " << solTemp.TE[i] << std::endl;
            }
            outStream << "y0 = \n";
            for(int i = 0;i < solTemp.Y[0].size();i++)
            {
                outStream << solTemp.Y[0][i] << std::endl;
            }
            outStream << "yend = \n";
            for(int i = 0;i < solTemp.Y[0].size();i++)
            {
                outStream << solTemp.Y[int(solTemp.Y.size() - 1)][i]
                          << std::endl;
            }
            outStream << "Tend = " << solTemp.T[int(solTemp.T.size()) - 1]
                      << std::endl;
            throw "Unexpected ode-event encountered during integration.";
        }
    }
    return overshoots;
}

//------------------------------------------------------------------------------
//Function to find the bounce for a single side (ie, don't attempt to find
//the other side, just do overshoot/undershoot on one side
template< class value_type , class time_type >
DLL_EXPORT void odeSolve_dS_one_sided(value_type false_vacuum ,
                     value_type true_vacuum ,
                     value_type barrier , potential< value_type >& V ,
                     time_type chimax_suggestion , int odeSolverToUse,
                     value_type RelTol, value_type AbsTol, value_type stepError,
                     value_type xi,
                     solution_grid< time_type , std::vector<value_type> ,
                     value_type >& solOut,value_type& DSout,
                     value_type precision,value_type h,value_type W0,
                     value_type& y_lower,value_type& y_upper,
                     std::ostream& outStream,bool track_scale_factor,
                     bool linear_BCs,value_type action_reltol,
                     value_type action_abstol,
                     bool graded_precision,
                     value_type a_scale,
                     value_type overshoot_lower,
                     value_type overshoot_upper,
                     bool userSuppliedBounds,
                     value_type LHS_bound_lower,
                     value_type LHS_bound_upper,
                     bool reverse_shooting,bool fixed_background)
{
    try
    {
//Instantiate the events function, state-type, dynamic BCs etc...
    typedef std::vector< value_type > dSSolType;//solution without action,
        //for speed.
    typedef std::vector< value_type > dSSolTypeAction;
    typedef ode_dS< dSSolType , time_type , value_type > odeRHS_dS;
    typedef eventsMethod0< dSSolType , time_type , value_type > ev_dS;
    typedef dynamicBCsAdSFlat_taylor< value_type , dSSolType , time_type >
            dynBC_dS;
    typedef ode_dS_action< dSSolTypeAction , time_type , value_type >
            odeRHS_dSAction;
    typedef eventsMethod0< dSSolTypeAction , time_type , value_type >
            ev_dSAction;
    typedef dynamicBCsAdSFlat_taylor< value_type , dSSolTypeAction ,
                                      time_type > dynBCAdSFlatAction;

    typedef eventsGravFriction<dSSolType,time_type,value_type> eventsGrav;
    typedef odeGravFriction< dSSolType , time_type , value_type > odeGrav;
    typedef dynamicBCsGravFriction_taylor< value_type , dSSolType ,
                                      time_type > dynBCGrav;
    typedef odeGravFrictionAction<dSSolType,time_type,value_type>
                                    odeGravAction;

    //Fixed background approximation ode functions:
    typedef ode_dS_fixed< dSSolType , time_type , value_type > odeRHS_dS_fixed;
    typedef ode_dS_fixed_action< dSSolType , time_type , value_type >
            odeRHS_dS_fixed_action;
    typedef events_dS_fixed< dSSolTypeAction , time_type , value_type >
            ev_dS_fixed;
    typedef dynamicBCs_dS_fixed< value_type , dSSolTypeAction ,
                                      time_type > dynBC_dS_fixed;




    odeRHS_dS odeToSolveTrack(V,xi, h, W0); //ode to solve (gravitational
                                                           //instanton equation)
    odeGrav odeToSolveNoTrack(V,xi, h, W0);
    odeRHS<std::vector<value_type>,time_type,value_type>* odeToSolve;


    //RHS_upper = true_vacuum;
    //Determine whether or not to use the user supplied bounds:
    const value_type _2p0 = value_type(2.0);
    const value_type _0p0 = value_type(0.0);
    const value_type _1p0 = value_type(1.0);
    const value_type _6p0 = value_type(6.0);
    const value_type _4p0 = value_type(4.0);
    const value_type _8p0 = value_type(8.0);
    const value_type _3p0 = value_type(3.0);
    const value_type _0p01 = value_type(0.01);
    const value_type _24p0 = value_type(24.0);
    const value_type _10p0 = value_type(10.0);
    const value_type PI = boost::math::constants::pi<value_type>();
    const value_type eps = std::numeric_limits<value_type>::epsilon();

    //Determine whether we have user supplied bounds, and if not, refine the
    //locations of the critical points and use those instead:
    value_type LHS_lower;
    value_type LHS_upper;
    if(userSuppliedBounds)
    {
        LHS_lower = LHS_bound_lower;
        LHS_upper = LHS_bound_upper;
    }
    else
    {
        LHS_lower = false_vacuum;
        LHS_upper = barrier;
    }

    //Remembers which side of the barrier we are one. 0 = false vacuum side,
    //1 = true vacuum side:
    int nSide;
    if((LHS_lower + LHS_upper)/_2p0 < barrier
        && (LHS_lower + LHS_upper)/_2p0 > false_vacuum)
    {
        nSide = 0;
    }
    else if((LHS_lower + LHS_upper)/_2p0 > barrier
        && (LHS_lower + LHS_upper)/_2p0 < true_vacuum)
    {
        nSide = 1;
    }
    else
    {
        throw "Error. Search must be confined to only one side of the barrier.";
    }
    std::cout << "\nnSide = " << nSide;



    /*
    value_type true_vacuum = tv;
    value_type barrier = bar;
    value_type LHS_lower = false_vacuum;
    value_type LHS_upper = bar;
    value_type RHS_lower = bar;
    value_type RHS_upper = tv;
    */
    std::cout << "\nfv = " << false_vacuum;
    std::cout << "\ntv = " << true_vacuum;

    //Specify the bounds used to identify overshoots (defaults to the false
    // and true vacuum)
    ev_dS eventsTrack(overshoot_lower,overshoot_upper);
    //events function.
    eventsGrav eventsNoTrack(overshoot_lower,overshoot_upper);
    events_function< std::vector<value_type> , time_type, value_type >* events;
    //ev_dS events(false_vacuum,true_vacuum); //events function.
        //Checks for events labelled IE:
        //IE = 1 -> a' crosses zero (terminal event, => overshoot)
        //IE = 2 -> y' crosses zero (terminal event, => undershoot)
        //IE = 3 -> a crosses zero (non-terminal event, for reference only)
        //IE = 4 -> y crosses upperBound (terminal event, indicates overshoot)
        //IE = 5 -> y crosses lowerBound (terminal event, indicates a->0)

    //Note - switch off termination for the event where we detect a'(x) = 0.
    //This will hinder the overshoot/undershoot method if we leave it on
    //because sometimes it isn't possible to tell whether a solution is an
    //overshoot or undershoot just by looking at the solution up until
    //a' = 0. Leave it on for the action version, however, because we only need
    //to compute half the solution from either side in that case anyway:
    eventsTrack.dynamicSwitch[0] = false;//Switches off termination for event 1
    eventsNoTrack.dynamicSwitch[0] = false;

    //Fixed background approximation:
    const value_type H0 = h*sqrt(W0/_3p0);
    odeRHS_dS_fixed_action odeToSolveFixedAction(V,xi, h, W0);
    ev_dS_fixed eventsActionFixed(false_vacuum,_2p0*true_vacuum,H0);
    eventsActionFixed.dynamicSwitch[0] = false;
    dynBC_dS_fixed bcAdjusterActionFixed(xi,h,W0,V,stepError,2);


    //Setup the functions which handle the boundary conditions at x = 0. This is
    //a co-ordinate singularity, so we need to perform a small step away from
    //it using a Taylor expansion in order to avoid numerical difficulties;

    //Two different versions, depending on whether we track a, a' separately,
    //('track' version) or just use an equation for a'/a ('no track' version).
    //Final argument determines whether we are using boundary conditions for the
    //action or not ('1' => no action, '2' => computing action)
    dynBC_dS bcAdjusterTrack(xi,h,W0,V,stepError,1);
    dynBCGrav bcAdjusterNoTrack(xi,h,W0,V,stepError,1);
    dynamicBCs_basicMethod< value_type, dSSolType, time_type >
        dynBCvac(xi,h,W0,V,stepError);//Applies to solutions computed without
            //action.
    dynamicBCs_actionMethod< value_type, dSSolTypeAction, time_type >
        dynBCvacAction(xi,h,W0,V,stepError);//applies to solutions computed
            //with action.
    //Setup pointers to the correct boundary conditions function, ode,
    //and events function, depending on whether
    //we use the 'track' or 'no track' version.
    dynamicBCs< value_type , std::vector<value_type> , time_type >* bcAdjuster;

    if(fixed_background)
    {
        std::cout << "\nUsing fixed background approximation.";
        odeToSolve = &odeToSolveFixedAction;
        bcAdjuster = &bcAdjusterActionFixed;
        events = &eventsActionFixed;
    }
    else if(track_scale_factor)
    {
        std::cout << "\nTracking full scale factor.";
        bcAdjuster = &bcAdjusterTrack;
        odeToSolve = &odeToSolveTrack;
        events = &eventsTrack;
    }
    else
    {
        std::cout << "\nUsing a'/a only.";
        bcAdjuster = &bcAdjusterNoTrack;
        odeToSolve = &odeToSolveNoTrack;
        events = &eventsNoTrack;
    }
    //Performs analytic first
    //step using a taylor series approximation.
        //stepError - determines size of initial step.
        //version (last argument) - specifies which version of the ode we are
            //using.
            //version = 1 -> normal ode, no action.
            //version = 2 -> ode with action





    //Search loop:
    int nMax = 200;//Terminate after this many iterations and generate a
        //a warning if we haven't found the bounce.
    int counter = 0;

    //Parameters related to checking whether we are trapped due to low
    //precision mis-identifying a bounce:
    value_type suspicion_thresh = value_type(0.001);//If we get a run with lower
        //probability than this, then we suspect that something has gone wrong
        //and switch to the maximum available precision.
    //Counters for measuring the a run of either type:
    int counter_undershoot_run = 0;
    int counter_overshoot_run = 0;
    //Resulting theshold value of N:
    value_type N_thresh = -log(suspicion_thresh)/log(_2p0);
    bool use_higher_precision = false;




    //Output some diagnostic parameters:
    //outStream << "Lower initial = " << lower << std::endl;
    //outStream << "Upper initial = " << upper << std::endl;

    outStream << "\nbarrier = " << barrier;
    outStream << "\nfv = " << false_vacuum;
    outStream << "\ntv = " << true_vacuum;

    outStream << "\nLHS_upper = " << LHS_upper;
    outStream << "\nLHS_lower = " << LHS_lower;

    //Difference between the bounds on the bounce (gives an estimate of how
    //cllose we are):
    value_type LHS_diff = abs(LHS_upper - LHS_lower);
    outStream << "LHS_diff = " << LHS_diff << std::endl;
    //Initial guess at bounce location:
    value_type LHS_guess = (LHS_lower + LHS_upper)/_2p0;


    //Tolerances to decide whether we are sufficiently close to the bounce
    //that we can stop looking:
    //value_type tol_LHS = precision*(abs(LHS_guess) + _1p0)/_2p0;
    //value_type tol_RHS = precision*(abs(RHS_guess) + _1p0)/_2p0;
    value_type tol_LHS = precision*abs(LHS_guess);
    outStream << "tol_LHS = " << tol_LHS << std::endl;

    //Maximum value of x to integrate up to:
    time_type chimax = chimax_suggestion;//Can add code here to adjust this
        //dynamically if we want, ie, take the user's input as a suggestion
        //only. May not be necessary, however.
    std::cout << "\nchimax = " << chimax;
    time_type initStepMax = time_type(0.01)/chimax;
    //Checks whether a(x) matches for the LHS and RHS solutions.
    bool a_match = false;

    //If we choose to steadily ramp up the precision as we get closer to the
    //bounce, use these initial values for the ode tolerances:
    value_type graded_reltol = RelTol;
    value_type graded_abstol = AbsTol;

    //struct to store the temporary results:
    solution_grid< time_type , dSSolType , value_type > LHS_solTemp;

    //Vector to store all the guesses we have made. This is useful for checking
    //whether errors have been made.
    std::vector<value_type> guess_vector;
    //Vector to store the reason that the loop terminated:
    std::vector<int> end_type;
    //Vector to store the terminating value of the solution (first component
                                                             //only).
    std::vector<value_type> end_value;
    std::vector<value_type> undershoot_list;
    std::vector<value_type> overshoot_list;

    //Add the initial bounds to the lists of overshoots and undershoots. These
    //are assumed to be accurate, so it is important that the user ensures this
    switch(nSide)
    {
    case 0://LHS of barrier (false vacuum side)
        if(reverse_shooting)
        {
            undershoot_list.push_back(LHS_lower);
            overshoot_list.push_back(LHS_upper);
        }
        else
        {
            overshoot_list.push_back(LHS_lower);
            undershoot_list.push_back(LHS_upper);
        }
        break;
    case 1://RHS of barrier (true vacuum side)
        if(reverse_shooting)
        {
            undershoot_list.push_back(LHS_upper);
            overshoot_list.push_back(LHS_lower);
        }
        else
        {
            overshoot_list.push_back(LHS_upper);
            undershoot_list.push_back(LHS_lower);
        }
        break;
    default:
        throw "Error - invalid side of potential.";
    }


    //Bounce search loop:
    do
    {
        //Integrate LHS solution (but don't bother if we already reached max
        //precision on it - just skip it and wait for the RHS to catch up).
        //solution_grid< time_type , dSSolType , value_type > LHS_solTemp;
        if(LHS_diff > tol_LHS)
        {
            LHS_solTemp.erase();
            //outStream << "Loop stage 1" <<std::endl;
            //Bisect to get initial guess:
            LHS_guess = (LHS_lower + LHS_upper)/_2p0;
            //tol_LHS = precision*(abs(LHS_guess) + _1p0)/_2p0;
            //tol_LHS = precision*abs(LHS_guess);
            tol_LHS = precision;
            outStream << "\nLHS_diff = " << LHS_diff;
            outStream << "\nLHS_y0 = " << LHS_guess;
            //outStream << "\nV'(LHS_y0) = " << V.d(LHS_guess);

            //Formulate initial conditions vector:
            dSSolType LHS_y0;
            if(!fixed_background)
            {
                value_type LHS_y0_List[] = { LHS_guess , _0p0 , _0p0 ,a_scale};
                LHS_y0.assign(LHS_y0_List,LHS_y0_List + 4);
            }
            else
            {
                value_type LHS_y0_List[] = { LHS_guess , _0p0 , _0p0 };
                LHS_y0.assign(LHS_y0_List,LHS_y0_List + 3);
            }

            //Step away to avoid co-ordinate singularity at chi = 0:
            //outStream << "Loop stage 3" <<std::endl;
            dSSolType LHS_y01step = LHS_y0;
            time_type LHS_epsilon;
            //if(false)
            if(abs(V.d(LHS_guess)) < stepError && linear_BCs)
            {
                LHS_epsilon = dynBCvac(LHS_y0,LHS_y01step);
            }
            else
            {
                LHS_epsilon = (*bcAdjuster)(LHS_y0,LHS_y01step);
            }
            for(int i = 0;i < LHS_y01step.size();i++)
            {
                outStream << "\nLHS_y01step[" << i << "] = " << LHS_y01step[i];
            }

            //Time span:
            time_type LHS_tspanArray[2] = { LHS_epsilon , chimax };
            std::vector< time_type > LHS_tspan;
            LHS_tspan.assign(LHS_tspanArray,LHS_tspanArray + 2);

            //Solve ode:
            outStream << "\nIntegrating LHS..." <<std::endl;
            if(graded_precision)
            {
                /*graded_reltol = std::min(value_type(1e-5)*
                                         abs(LHS_diff/LHS_guess),RelTol);*/
                graded_reltol = std::min(abs(LHS_diff),value_type(1e-10));
                graded_abstol = AbsTol;
            }
            else
            {
                graded_reltol = RelTol;
                graded_abstol = AbsTol;
            }
            //Override and use higher precision if we suspect an error in
            //solution classification:
            if(use_higher_precision)
            {
                graded_reltol = action_reltol;
                graded_abstol = action_abstol;
            }
            std::cout << "\ngraded_reltol = " << graded_reltol;
            int LHS_nSuccess = ode_solve<time_type,dSSolType,value_type>
                (*odeToSolve,LHS_y01step,LHS_tspan,*events,odeSolverToUse,
                 LHS_solTemp,graded_reltol,graded_abstol,initStepMax);
            if(LHS_nSuccess != 1)
            {
                outStream << "nSuccess = " << LHS_nSuccess
                          << " when solving for AdS-Flat instantons."
                          << std::endl;
                solOut.Y = LHS_solTemp.Y;
                solOut.T = LHS_solTemp.T;
                solOut.YE = LHS_solTemp.YE;
                solOut.TE = LHS_solTemp.TE;
                solOut.IE = LHS_solTemp.IE;
                throw "Integration failure.";
            }

            //Now check whether we found an overshoot or an undershoot:
            bool LHS_overshoots = check_overshoot(LHS_solTemp,LHS_upper,
                                                  LHS_lower,
                                                  LHS_guess,outStream,nSide,
                                                  reverse_shooting);
            std::cout << "\n" << (LHS_overshoots ? "overshoot " : "undershoot ")
                      << "\nLHS_lower = " << LHS_lower
                      << "\nLHS_upper = " << LHS_upper
                      << "\nLHS_guess = " << LHS_guess
                      << "\nIE = (";
            for(int i = 0;i < LHS_solTemp.IE.size();i++)
            {
                if(i == 0)
                {
                    std::cout << LHS_solTemp.IE[i];
                }
                else
                {
                    std::cout << "," << LHS_solTemp.IE[i];
                }
            }
            std::cout << ").\n";
            //Record solution's initial value:
            guess_vector.push_back(LHS_guess);
            //Record why the solution terminated:
            int nEndValue = LHS_solTemp.IE[LHS_solTemp.IE.size() - 1];
            end_type.push_back(nEndValue);
            //And its final value:
            end_value.push_back(LHS_solTemp.Y[LHS_solTemp.Y.size() - 1][0]);
            //Store a list of overshoot and undershoots:
            if(nEndValue == 2)
            {
                undershoot_list.push_back(LHS_guess);
            }
            else
            {
                overshoot_list.push_back(LHS_guess);
            }

            //Check that we don't get too many of one category in a row, as this
            //may indicate that something has gone wrong and the algorithm is
            //trapped:
            if(LHS_overshoots)
            {
                counter_overshoot_run ++;
                counter_undershoot_run = 0;
            }
            else
            {
                counter_undershoot_run ++;
                counter_overshoot_run = 0;
            }

            if(value_type(counter_undershoot_run) > N_thresh ||
               value_type(counter_overshoot_run) > N_thresh )
            {
                if(!use_higher_precision)
                {
                    use_higher_precision = true;
                    //Go back through the list to check the undershoots - one
                    //of them might be wrong!
                    std::cout << "\nMiss-classification suspected."
                        << " Check at higher precision: ";
                    int current_end_type = LHS_solTemp.IE[LHS_solTemp.IE.size()
                                                          - 1];
                    for(int i = overshoot_list.size() - 1;i >= 0;i--)
                    {
                        LHS_guess = overshoot_list[i];
                        dSSolType y0;
                        if(!fixed_background)
                        {
                            value_type y0_List[] = {LHS_guess,_0p0
                                                      ,_0p0,a_scale};
                            y0.assign(y0_List,y0_List + 4);
                        }
                        else
                        {
                            value_type y0_List[] = {LHS_guess,_0p0,_0p0};
                            y0.assign(y0_List,y0_List + 3);
                        }


                        dSSolType y01step = y0;
                        time_type eps;

                        std::cout.precision(50);
                        std::cout << "\nY_{lhs} = ";
                        for(int i = 0;i < y0.size();i++)
                        {
                            std::cout << "\nY[" << i << "] = "
                                      << y0[i];
                        }

                        //Adjust boundary conditions to avoid co-ordinate
                        //singularities:
                        //if(false)
                        if(abs(V.d(LHS_guess)) < stepError && linear_BCs)
                        {
                            eps = dynBCvac(y0,y01step);
                        }
                        else
                        {
                            eps = (*bcAdjuster)(y0,y01step);
                        }

                        //Clear solution grids for the final result:
                        LHS_solTemp.erase();


                        //Time span to integrate over:
                        time_type tspanArray[2] = { eps , chimax };
                        std::vector< time_type > tspan;
                        tspan.assign(tspanArray,tspanArray + 2);
                        initStepMax = time_type(0.01)/chimax;

                        //Since we are only computing the bounce in one
                        //direction, switch
                        //off termination for a' = 0 event;
                        eventsTrack.dynamicSwitch[0] = false;
                        eventsNoTrack.dynamicSwitch[0] = false;



                        //LHS solution:
                        std::cout << "\nComputing bounce solution.";
                        std::cout << "\naction_reltol = " << action_reltol;
                        std::cout << "\naction_abstol = " << action_abstol;
                        int nSuccess = ode_solve<time_type,dSSolType,value_type>
                        (*odeToSolve,y01step,tspan,*events,odeSolverToUse,
                         LHS_solTemp,action_reltol,action_abstol,initStepMax);
                        if(nSuccess != 1)
                        {
                            outStream << "nSuccess = " << LHS_nSuccess
                                      << " when solving for AdS-Flat"
                                      << "instantons." << std::endl;
                            throw "Integration failure.";
                        }
                        std::cout << "\nComputed " << LHS_solTemp.T.size()
                                  << " points on the LHS.";
                        std::cout << "\nEvents = (";
                        for(int i = 0;i < LHS_solTemp.IE.size();i++)
                        {
                            if(i == 0)
                            {
                                std::cout << LHS_solTemp.IE[i];
                            }
                            else
                            {
                                std::cout << "," << LHS_solTemp.IE[i];
                            }
                        }
                        std::cout << ").\n";

                        //Check whether this worked:
                        current_end_type = LHS_solTemp.IE[LHS_solTemp.IE.size()
                                                          - 1];
                        if(current_end_type == 4 || current_end_type == 5)
                        {
                            //This code will re-assign the relevant bounds,
                            //probably widening them if we mis-identified
                            //anything.
                            LHS_overshoots = check_overshoot(LHS_solTemp,
                                                                  LHS_upper,
                                                  LHS_lower,
                                                  LHS_guess,outStream,nSide,
                                                  reverse_shooting);
                            break;
                        }
                    }
                    std::cout << "\nCompleted checks on overshoots.";
                    //Once here, we either found a reliable overshoot, or
                    //we got back to the initial bounds.
                    //Go back through the list to check the undershoots - one
                    //of them might be wrong!
                    current_end_type = LHS_solTemp.IE[LHS_solTemp.IE.size()
                                                          - 1];
                    for(int i = undershoot_list.size() - 1;i >= 0;i--)
                    {
                        LHS_guess = undershoot_list[i];
                        dSSolType y0;
                        if(!fixed_background)
                        {
                            value_type y0_List[] ={LHS_guess,_0p0,_0p0,a_scale};
                            y0.assign(y0_List,y0_List + 4);
                        }
                        else
                        {
                            value_type y0_List[] ={LHS_guess,_0p0,_0p0};
                            y0.assign(y0_List,y0_List + 3);
                        }
                        dSSolType y01step = y0;
                        time_type eps;

                        std::cout.precision(50);
                        std::cout << "\nY_{lhs} = ";
                        for(int i = 0;i < y0.size();i++)
                        {
                            std::cout << "\nY[" << i << "] = "
                                      << y0[i];
                        }

                        //Adjust boundary conditions to avoid co-ordinate
                        //singularities:
                        //if(false)
                        if(abs(V.d(LHS_guess)) < stepError && linear_BCs)
                        {
                            eps = dynBCvac(y0,y01step);
                        }
                        else
                        {
                            eps = (*bcAdjuster)(y0,y01step);
                        }

                        //Clear solution grids for the final result:
                        LHS_solTemp.erase();


                        //Time span to integrate over:
                        time_type tspanArray[2] = { eps , chimax };
                        std::vector< time_type > tspan;
                        tspan.assign(tspanArray,tspanArray + 2);
                        initStepMax = time_type(0.01)/chimax;

                        //Since we are only computing the bounce in one
                        //direction, switch
                        //off termination for a' = 0 event;
                        eventsTrack.dynamicSwitch[0] = false;
                        eventsNoTrack.dynamicSwitch[0] = false;



                        //LHS solution:
                        std::cout << "\nComputing bounce solution.";
                        std::cout << "\naction_reltol = " << action_reltol;
                        std::cout << "\naction_abstol = " << action_abstol;
                        int nSuccess = ode_solve<time_type,dSSolType,value_type>
                        (*odeToSolve,y01step,tspan,*events,odeSolverToUse,
                         LHS_solTemp,action_reltol,action_abstol,initStepMax);
                        if(nSuccess != 1)
                        {
                            outStream << "nSuccess = " << LHS_nSuccess
                                      << " when solving for AdS-Flat"
                                      << "instantons." << std::endl;
                            throw "Integration failure.";
                        }
                        std::cout << "\nComputed " << LHS_solTemp.T.size()
                                  << " points on the LHS.";
                        std::cout << "\nEvents = (";
                        for(int i = 0;i < LHS_solTemp.IE.size();i++)
                        {
                            if(i == 0)
                            {
                                std::cout << LHS_solTemp.IE[i];
                            }
                            else
                            {
                                std::cout << "," << LHS_solTemp.IE[i];
                            }
                        }
                        std::cout << ").\n";

                        //Check whether this worked:
                        current_end_type = LHS_solTemp.IE[LHS_solTemp.IE.size()
                                                          - 1];
                        if(current_end_type == 2)
                        {
                            //This code will re-assign the relevant bounds,
                            //probably widening them if we mis-identified
                            //anything.
                            LHS_overshoots = check_overshoot(LHS_solTemp,
                                                                  LHS_upper,
                                                  LHS_lower,
                                                  LHS_guess,outStream,nSide,
                                                  reverse_shooting);
                            //break;
                            break;
                        }
                    }
                    std::cout << "\nCompleted checks on undershoots.";
                    //Once here, we either found a reliable undershoot, or
                    //we got back to the initial bounds.
                    std::cout << "\nNew bounds: ";
                    std::cout << "\nLower = " << LHS_lower;
                    std::cout << "\nUpper = " << LHS_upper;
                    //Reset counters, otherwise we get erroneous complaints
                    //about runs of too many overshoot/undershoots:
                    counter_undershoot_run = 0;
                    counter_overshoot_run = 0;
                }
                else
                {
                    //Exit loop - we can't proceed any further:
                    std::cout << "\nSuspected miss-classification,but "
                              << "reached maximum accuracy. Results may "
                              << "be inaccurate.";
                    break;
                }
            }
        }


        counter++;
        LHS_diff = abs(LHS_upper - LHS_lower);
        std::cout << "\nLHS_diff = " << LHS_diff;
        std::cout << "\nLHS_lower = " << LHS_lower;
        std::cout << "\nLHS_upper = " << LHS_upper;
        if(counter > nMax)
        {
            break;
        }

        //Clear solution grids for the next pass, unless we are done with
        //that half:
        if(LHS_diff > tol_LHS)
        {
            LHS_solTemp.erase();
        }
    }
    while(LHS_diff > tol_LHS);


    //We found the bounce. Are we within precision?
    if(LHS_diff > tol_LHS)
    {
        outStream << "Loop exited before finding bounce." << std::endl;
        //reportToCaller("Loop exited before finding bounce.\n");
    }

    //Now compute the action. This is typically done at much higher precision:
    outStream << "\nSetup action version.";
    odeRHS_dSAction odeToSolveActionTrack(V,xi,h,W0);
    //outStream << "\n1";
    ev_dSAction eventsActionTrack(false_vacuum,true_vacuum);
    //outStream << "\n2";
    dynBCAdSFlatAction bcAdjusterActionTrack(xi,h,W0,V,stepError,2);
    //outStream << "\n3";
    odeGravAction odeToSolveActionNoTrack(V,xi,h,W0);
    //outStream << "\n4";
    eventsGrav eventsActionNoTrack(false_vacuum,true_vacuum);
    //outStream << "\n5";
    dynBCGrav bcAdjusterActionNoTrack(xi,h,W0,V,stepError,2);
    //outStream << "\n6";

    dynamicBCs< value_type , std::vector<value_type> , time_type >*
        bcAdjusterAction;
    //outStream << "\n7";
    events_function< std::vector<value_type> , time_type, value_type >*
        eventsAction;
    //outStream << "\n8";
    odeRHS<std::vector<value_type>,time_type,value_type>* odeToSolveAction;
    //outStream << "\n9";
    if(fixed_background)
    {
        //std::cout << "\nUsing fixed background approximation.";
        odeToSolveAction = &odeToSolveFixedAction;
        bcAdjusterAction = &bcAdjusterActionFixed;
        eventsAction = &eventsActionFixed;
    }
    else if(track_scale_factor)
    {
        odeToSolveAction = &odeToSolveActionTrack;
        bcAdjusterAction = &bcAdjusterActionTrack;
        eventsAction = &eventsActionTrack;
    }
    else
    {
        odeToSolveAction = &odeToSolveActionNoTrack;
        bcAdjusterAction = &bcAdjusterActionNoTrack;
        eventsAction = &eventsActionNoTrack;
    }
    //outStream << "\n10";

    //Bounce boundary conditions:
    dSSolTypeAction LHS_y0Action;
    //LHS_guess = (LHS_upper + LHS_lower)/_2p0;
    //We prefer undershoots to overshoots, as these are
    //generally terminated before they diverge, so lie closer to the bounce
    //solution (in general):
    if(nSide == 0)
    {
        if(reverse_shooting)
        {
            LHS_guess = LHS_lower;
        }
        else
        {
            LHS_guess = LHS_upper;
        }
    }
    else
    {
        if(reverse_shooting)
        {
            LHS_guess = LHS_upper;
        }
        else
        {
            LHS_guess = LHS_lower;
        }
    }
    if(!fixed_background)
    {
        value_type LHS_y0Action_List[] = {LHS_guess,_0p0,_0p0,a_scale,_0p0};
        LHS_y0Action.assign(LHS_y0Action_List,LHS_y0Action_List + 5);
    }
    else
    {
        value_type LHS_y0Action_List[] = {LHS_guess,_0p0,_0p0};
        LHS_y0Action.assign(LHS_y0Action_List,LHS_y0Action_List + 3);
    }

    dSSolTypeAction LHS_y0Action1step = LHS_y0Action;
    time_type LHS_epsAction;

    std::cout.precision(50);
    std::cout << "\nY_{lhs} = ";
    for(int i = 0;i < LHS_y0Action.size();i++)
    {
        std::cout << "\nY[" << i << "] = " << LHS_y0Action[i];
    }

    //Adjust boundary conditions to avoid co-ordinate singularities:
    //if(false)
    if(abs(V.d(LHS_guess)) < stepError && linear_BCs)
    {
        LHS_epsAction = dynBCvacAction(LHS_y0Action,LHS_y0Action1step);
    }
    else
    {
        LHS_epsAction = (*bcAdjusterAction)(LHS_y0Action,LHS_y0Action1step);
    }

    //Clear solution grids for the final result:
    LHS_solTemp.erase();


    //Time span to integrate over:
    time_type LHS_tspanArray[2] = { LHS_epsAction , chimax };
    std::vector< time_type > LHS_tspan;
    LHS_tspan.assign(LHS_tspanArray,LHS_tspanArray + 2);
    initStepMax = time_type(0.01)/chimax;

    //Since we are only computing the bounce in one direction, switch
    //off termination for a' = 0 event;
    eventsActionTrack.dynamicSwitch[0] = false;
    eventsActionNoTrack.dynamicSwitch[0] = false;



    //LHS solution:
    std::cout << "\nComputing bounce solution.";
    std::cout << "\naction_reltol = " << action_reltol;
    std::cout << "\naction_abstol = " << action_abstol;
    int LHS_nSuccess = ode_solve<time_type,dSSolTypeAction,value_type>
    (*odeToSolveAction,LHS_y0Action1step,LHS_tspan,*eventsAction,odeSolverToUse,
     LHS_solTemp,action_reltol,action_abstol,initStepMax);
    if(LHS_nSuccess != 1)
    {
        outStream << "nSuccess = " << LHS_nSuccess
                  << " when solving for AdS-Flat"
                  << "instantons." << std::endl;
        throw "Integration failure.";
    }
    std::cout << "\nComputed " << LHS_solTemp.T.size() << " points on the LHS.";
    std::cout << "\nEvents = (";
    for(int i = 0;i < LHS_solTemp.IE.size();i++)
    {
        if(i == 0)
        {
            std::cout << LHS_solTemp.IE[i];
        }
        else
        {
            std::cout << "," << LHS_solTemp.IE[i];
        }
    }
    std::cout << ").\n";

    //Output a list of all the solutions we found:
    std::cout.precision(50);
    std::cout << "\nSolutions considered in search: ";
    for(int i = 0;i < guess_vector.size();i++)
    {
        std::cout << "\ny0 = " << guess_vector[i];
        std::cout << "\n"
            << (end_type[i] == 2 ? "undershoot: IE = " : "overshoot: IE = ")
            << end_type[i];
        std::cout << "\ny_end = " << end_value[i] << std::endl;
    }

    //Check whether the solution is an undershoot, as we require:
    if(LHS_solTemp.IE[LHS_solTemp.IE.size() - 1] != 2)
    {
        //Algorithm got trapped, probably due to a miss-classification
        //somewhere along the line. Go back through the solutions we computed
        //and attempt to find the correct location.
        std::cout << "\n\nSolution is not an undershoot. Attempt to correct: ";
        int current_end_type = LHS_solTemp.IE[LHS_solTemp.IE.size() - 1];
        for(int i = undershoot_list.size() - 1;i >= 0;i--)
        {
            dSSolTypeAction y0Action;
            if(!fixed_background)
            {
                value_type y0Action_List[] = {undershoot_list[i],_0p0,_0p0,
                                          a_scale,_0p0};
                y0Action.assign(y0Action_List,y0Action_List + 5);
            }
            else
            {
                value_type y0Action_List[] = {undershoot_list[i],_0p0,_0p0};
                y0Action.assign(y0Action_List,y0Action_List + 3);
            }
            dSSolTypeAction y0Action1step = y0Action;
            time_type epsAction;

            std::cout.precision(50);
            std::cout << "\nY_{lhs} = ";
            for(int i = 0;i < y0Action.size();i++)
            {
                std::cout << "\nY[" << i << "] = " << y0Action[i];
            }

            //Adjust boundary conditions to avoid co-ordinate singularities:
            //if(false)
            if(abs(V.d(LHS_guess)) < stepError && linear_BCs)
            {
                epsAction = dynBCvacAction(y0Action,y0Action1step);
            }
            else
            {
                epsAction = (*bcAdjusterAction)(y0Action,y0Action1step);
            }

            //Clear solution grids for the final result:
            LHS_solTemp.erase();


            //Time span to integrate over:
            time_type tspanArray[2] = { epsAction , chimax };
            std::vector< time_type > tspan;
            tspan.assign(tspanArray,tspanArray + 2);
            initStepMax = time_type(0.01)/chimax;

            //Since we are only computing the bounce in one direction, switch
            //off termination for a' = 0 event;
            eventsActionTrack.dynamicSwitch[0] = false;
            eventsActionNoTrack.dynamicSwitch[0] = false;



            //LHS solution:
            std::cout << "\nComputing bounce solution.";
            std::cout << "\naction_reltol = " << action_reltol;
            std::cout << "\naction_abstol = " << action_abstol;
            int nSuccess = ode_solve<time_type,dSSolTypeAction,value_type>
            (*odeToSolveAction,y0Action1step,tspan,*eventsAction,odeSolverToUse,
             LHS_solTemp,action_reltol,action_abstol,initStepMax);
            if(nSuccess != 1)
            {
                outStream << "nSuccess = " << LHS_nSuccess
                          << " when solving for AdS-Flat"
                          << "instantons." << std::endl;
                throw "Integration failure.";
            }
            std::cout << "\nComputed " << LHS_solTemp.T.size()
                      << " points on the LHS.";
            std::cout << "\nEvents = (";
            for(int i = 0;i < LHS_solTemp.IE.size();i++)
            {
                if(i == 0)
                {
                    std::cout << LHS_solTemp.IE[i];
                }
                else
                {
                    std::cout << "," << LHS_solTemp.IE[i];
                }
            }
            std::cout << ").\n";

            //Check whether this worked:
            current_end_type = LHS_solTemp.IE[LHS_solTemp.IE.size() - 1];
            if(current_end_type == 2)
            {
                break;
            }
        }
    }
    //Output the bounds so that we can compute a higher solution is we want:
    std::cout << "\nLower bound = " << LHS_lower;
    std::cout << "\nUpper bound = " << LHS_upper;

    outStream << "\nSolution Computed.";

    //Output the solution:
    //solOut = LHS_solTemp;
    for(int i = 0;i < LHS_solTemp.T.size();i++)
    {
        solOut.T.push_back(LHS_solTemp.T[i]);
        std::vector<value_type> value_i (5,_0p0);
        value_i[0] = LHS_solTemp.Y[i][0];
        value_i[1] = LHS_solTemp.Y[i][1];
        value_i[2] = sin(H0*value_type(LHS_solTemp.T[i]))/H0;
        value_i[3] = cos(H0*value_type(LHS_solTemp.T[i]));
        value_i[4] = LHS_solTemp.Y[i][2];
        solOut.Y.push_back(value_i);
    }
    for(int i = 0;i < LHS_solTemp.TE.size();i++)
    {
        std::vector<value_type> value_i (5,_0p0);
        value_i[0] = LHS_solTemp.YE[i][0];
        value_i[1] = LHS_solTemp.YE[i][1];
        value_i[2] = sin(H0*value_type(LHS_solTemp.TE[i]))/H0;
        value_i[3] = cos(H0*value_type(LHS_solTemp.TE[i]));
        value_i[4] = LHS_solTemp.YE[i][2];
        solOut.YE.push_back(value_i);
        solOut.TE.push_back(LHS_solTemp.TE[i]);
        solOut.IE.push_back(LHS_solTemp.IE[i]);
    }
    DSout = !fixed_background ? solOut.Y[int(solOut.T.size()) - 1][4] :
                                solOut.Y[int(solOut.T.size()) - 1][2];

    //Record the upper and lower bounds:
    y_lower = LHS_lower;
    y_upper = LHS_upper;


    outStream << "Done." << std::endl;
    //Done.
    }
    catch(char* c)
    {
        outStream << c;
    }
    catch(...)
    {
        outStream << "An unknown error occurred.";
    }
}

//------------------------------------------------------------------------------
//Overloaded versions of ode_solve_dS using default parameters:
template< class value_type , class time_type >
DLL_EXPORT void odeSolve_dS(value_type false_vacuum , value_type true_vacuum ,
                     value_type barrier , potential< value_type >& V ,
                     time_type chimax_suggestion , int odeSolverToUse,
                     value_type RelTol, value_type AbsTol, value_type stepError,
                     value_type xi,
                     solution_grid< time_type , std::vector<value_type> ,
                     value_type >& solOut,
                     solution_grid< time_type , std::vector< value_type > ,
                     value_type >& solOutLHS,
                     solution_grid< time_type , std::vector< value_type > ,
                     value_type >& solOutRHS,
                     value_type& DSout,
                     value_type precision,value_type h,value_type W0,
                     std::vector<value_type>& other_return_values,
                     std::ostream& outStream,bool track_scale_factor,
                     bool linear_BCs)
{
    //Calls default parameters:
    /*
    action_rel_tol = RelTol;
    action_abs_tol = AbsTol;
    graded_precision = false;
    a_scale = 1.0;
    overshoot_lower = false_vacuum
    overshoot_upper = true_vacuum
    userSuppliedBounds = false
    LHS_bound_lower = false_vacuum
    LHS_bound_upper = barrier
    RHS_bound_lower = barrier
    RHS_bound_upper = true_vacuum
    reverse_shooting = false
    use_analytic_a = false
    */
    odeSolve_dS(false_vacuum , true_vacuum ,  barrier , V , chimax_suggestion ,
                odeSolverToUse, RelTol, AbsTol, stepError,xi,solOut,solOutLHS,
                solOutRHS,DSout,precision,h,W0,other_return_values,outStream,
                track_scale_factor,linear_BCs = false,RelTol,AbsTol,false,
                value_type(1.0),false_vacuum,true_vacuum,false,false_vacuum,
                barrier,barrier,true_vacuum,false,false,false);
}
//------------------------------------------------------------------------------
template< class value_type , class time_type >
DLL_EXPORT void odeSolve_dS(value_type false_vacuum , value_type true_vacuum ,
                     value_type barrier , potential< value_type >& V ,
                     time_type chimax_suggestion , int odeSolverToUse,
                     value_type RelTol, value_type AbsTol, value_type stepError,
                     value_type xi,
                     solution_grid< time_type , std::vector<value_type> ,
                     value_type >& solOut,
                     solution_grid< time_type , std::vector< value_type > ,
                     value_type >& solOutLHS,
                     solution_grid< time_type , std::vector< value_type > ,
                     value_type >& solOutRHS,
                     value_type& DSout,
                     value_type precision,value_type h,value_type W0,
                     std::vector<value_type>& other_return_values,
                     std::ostream& outStream,bool track_scale_factor,
                     bool linear_BCs,value_type action_reltol,
                     value_type action_abstol,bool graded_precision,
                     value_type a_scale)
{
    //Calls default parameters:
    /*
    overshoot_lower = false_vacuum
    overshoot_upper = true_vacuum
    userSuppliedBounds = false
    LHS_bound_lower = false_vacuum
    LHS_bound_upper = barrier
    RHS_bound_lower = barrier
    RHS_bound_upper = true_vacuum
    reverse_shooting = false
    use_analytic_a = false
    */
    odeSolve_dS(false_vacuum , true_vacuum ,  barrier , V , chimax_suggestion ,
                odeSolverToUse, RelTol, AbsTol, stepError,xi,solOut,solOutLHS,
                solOutRHS,DSout,precision,h,W0,other_return_values,outStream,
                track_scale_factor,linear_BCs = false,action_reltol,
                action_abstol,graded_precision,a_scale,false_vacuum,
                true_vacuum,false,false_vacuum,
                barrier,barrier,true_vacuum,false,false,false);
}
//------------------------------------------------------------------------------
template< class value_type , class time_type >
DLL_EXPORT void odeSolve_dS(value_type false_vacuum , value_type true_vacuum ,
                     value_type barrier , potential< value_type >& V ,
                     time_type chimax_suggestion , int odeSolverToUse,
                     value_type RelTol, value_type AbsTol, value_type stepError,
                     value_type xi,
                     solution_grid< time_type , std::vector<value_type> ,
                     value_type >& solOut,
                     solution_grid< time_type , std::vector< value_type > ,
                     value_type >& solOutLHS,
                     solution_grid< time_type , std::vector< value_type > ,
                     value_type >& solOutRHS,
                     value_type& DSout,
                     value_type precision,value_type h,value_type W0,
                     std::vector<value_type>& other_return_values,
                     std::ostream& outStream,bool track_scale_factor,
                     bool linear_BCs,value_type action_reltol,
                     value_type action_abstol,bool graded_precision,
                     value_type a_scale,value_type overshoot_lower,
                     value_type overshoot_upper,bool userSuppliedBounds)
{
    //Calls default parameters:
    /*
    LHS_bound_lower = false_vacuum
    LHS_bound_upper = barrier
    RHS_bound_lower = barrier
    RHS_bound_upper = true_vacuum
    reverse_shooting = false
    use_analytic_a = false
    */
    odeSolve_dS(false_vacuum , true_vacuum ,  barrier , V , chimax_suggestion ,
                odeSolverToUse, RelTol, AbsTol, stepError,xi,solOut,solOutLHS,
                solOutRHS,DSout,precision,h,W0,other_return_values,outStream,
                track_scale_factor,linear_BCs = false,action_reltol,
                action_abstol,graded_precision,a_scale,overshoot_lower,
                overshoot_upper,userSuppliedBounds,false_vacuum,
                barrier,barrier,true_vacuum,false,false,false);
}
//------------------------------------------------------------------------------
template< class value_type , class time_type >
DLL_EXPORT void odeSolve_dS_one_sided(value_type false_vacuum ,
                     value_type true_vacuum ,
                     value_type barrier , potential< value_type >& V ,
                     time_type chimax_suggestion , int odeSolverToUse,
                     value_type RelTol, value_type AbsTol, value_type stepError,
                     value_type xi,
                     solution_grid< time_type , std::vector<value_type> ,
                     value_type >& solOut,value_type& DSout,
                     value_type precision,value_type h,value_type W0,
                     value_type& y_lower,value_type& y_upper,
                     std::ostream& outStream,bool track_scale_factor,
                     bool linear_BCs)
{
    /*
    Default parameters used:
    action_reltol = RelTol
    action_abstol = AbsTol
    graded_precision = false
    a_scale = value_type(1.0)
    overshoot_lower = false_vacuum
    overshoot_upper = true_vacuum
    userSuppliedBounds = false
    LHS_bound_lower = false_vacuum
    LHS_bound_upper = barrier
    reverse_shooting = false
    */
    odeSolve_dS_one_sided(false_vacuum , true_vacuum ,barrier , V ,
                          chimax_suggestion , odeSolverToUse,
                          RelTol, AbsTol, stepError,
                          xi,solOut,DSout,precision,h,W0,
                          y_lower,y_upper,outStream,track_scale_factor,
                          linear_BCs,RelTol,AbsTol,false,
                          value_type(1.0),false_vacuum,true_vacuum,false,
                          false_vacuum,barrier,false,false);
}
//------------------------------------------------------------------------------
template< class value_type , class time_type >
DLL_EXPORT void odeSolve_dS_one_sided(value_type false_vacuum ,
                     value_type true_vacuum ,
                     value_type barrier , potential< value_type >& V ,
                     time_type chimax_suggestion , int odeSolverToUse,
                     value_type RelTol, value_type AbsTol, value_type stepError,
                     value_type xi,
                     solution_grid< time_type , std::vector<value_type> ,
                     value_type >& solOut,value_type& DSout,
                     value_type precision,value_type h,value_type W0,
                     value_type& y_lower,value_type& y_upper,
                     std::ostream& outStream,bool track_scale_factor,
                     bool linear_BCs,value_type action_reltol,
                     value_type action_abstol,bool graded_precision,
                     value_type a_scale)
{
    /*
    Default parameters used:
    overshoot_lower = false_vacuum
    overshoot_upper = true_vacuum
    userSuppliedBounds = false
    LHS_bound_lower = false_vacuum
    LHS_bound_upper = barrier
    reverse_shooting = false
    */
    odeSolve_dS_one_sided(false_vacuum , true_vacuum ,barrier , V ,
                          chimax_suggestion , odeSolverToUse,
                          RelTol, AbsTol, stepError,
                          xi,solOut,DSout,precision,h,W0,
                          y_lower,y_upper,outStream,track_scale_factor,
                          linear_BCs,action_reltol,action_abstol,
                          graded_precision,a_scale,false_vacuum,true_vacuum,
                          false,false_vacuum,barrier,false,false);
}
//------------------------------------------------------------------------------
template< class value_type , class time_type >
DLL_EXPORT void odeSolve_dS_one_sided(value_type false_vacuum ,
                     value_type true_vacuum ,
                     value_type barrier , potential< value_type >& V ,
                     time_type chimax_suggestion , int odeSolverToUse,
                     value_type RelTol, value_type AbsTol, value_type stepError,
                     value_type xi,
                     solution_grid< time_type , std::vector<value_type> ,
                     value_type >& solOut,value_type& DSout,
                     value_type precision,value_type h,value_type W0,
                     value_type& y_lower,value_type& y_upper,
                     std::ostream& outStream,bool track_scale_factor,
                     bool linear_BCs,value_type action_reltol,
                     value_type action_abstol,
                     bool graded_precision,
                     value_type a_scale,
                     value_type overshoot_lower,
                     value_type overshoot_upper,
                     bool userSuppliedBounds)
{
    /*
    Default parameters used:
    LHS_bound_lower = false_vacuum
    LHS_bound_upper = barrier
    reverse_shooting = false
    */
    odeSolve_dS_one_sided(false_vacuum , true_vacuum ,barrier , V ,
                          chimax_suggestion , odeSolverToUse,
                          RelTol, AbsTol, stepError,
                          xi,solOut,DSout,precision,h,W0,
                          y_lower,y_upper,outStream,track_scale_factor,
                          linear_BCs,action_reltol,action_abstol,
                          graded_precision,a_scale,overshoot_lower,
                          overshoot_upper,userSuppliedBounds,false_vacuum,
                          barrier,false,false);
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//Explicit instantiation and export of this function:
template DLL_EXPORT void odeSolve_dS
    < multi , multi >
    (multi , multi , multi , potential< multi >& , multi , int,multi, multi ,
     multi, multi, solution_grid< multi ,  std::vector<multi> , multi >&,
     solution_grid< multi ,  std::vector<multi> , multi >&,
     solution_grid< multi ,  std::vector<multi> , multi >&,
     multi&, multi, multi, multi,std::vector<multi>&,std::ostream&,bool,bool,
     multi,multi,bool,multi,multi,multi,bool,multi,multi,multi,multi,bool,bool,
     bool);
//------------------------------------------------------------------------------
template DLL_EXPORT void odeSolve_dS
    < multi , multi >
    (multi , multi , multi , potential< multi >& , multi , int,multi, multi ,
     multi, multi, solution_grid< multi ,  std::vector<multi> , multi >&,
     solution_grid< multi ,  std::vector<multi> , multi >&,
     solution_grid< multi ,  std::vector<multi> , multi >&,
     multi&, multi, multi, multi,std::vector<multi>&,std::ostream&,bool,bool);
//------------------------------------------------------------------------------
template DLL_EXPORT void odeSolve_dS
    < multi , multi >
    (multi , multi , multi , potential< multi >& , multi , int,multi, multi ,
     multi, multi, solution_grid< multi ,  std::vector<multi> , multi >&,
     solution_grid< multi ,  std::vector<multi> , multi >&,
     solution_grid< multi ,  std::vector<multi> , multi >&,
     multi&, multi, multi, multi,std::vector<multi>&,std::ostream&,bool,bool,
     multi,multi,bool,multi);
//------------------------------------------------------------------------------
template DLL_EXPORT void odeSolve_dS
    < multi , multi >
    (multi , multi , multi , potential< multi >& , multi , int,multi, multi ,
     multi, multi, solution_grid< multi ,  std::vector<multi> , multi >&,
     solution_grid< multi ,  std::vector<multi> , multi >&,
     solution_grid< multi ,  std::vector<multi> , multi >&,
     multi&, multi, multi, multi,std::vector<multi>&,std::ostream&,bool,bool,
     multi,multi,bool,multi,multi,multi,bool);
//------------------------------------------------------------------------------
template DLL_EXPORT void odeSolve_dSsingle
    < multi , multi >
    (multi , multi , multi , multi , potential< multi >& , multi , int, multi,
     multi , multi, multi, solution_grid< multi , std::vector<multi> , multi >&,
     multi&, multi, multi, multi, std::ostream&,bool,bool,bool);
//------------------------------------------------------------------------------
template DLL_EXPORT void odeSolve_dS_one_sided
    < multi , multi >
    (multi , multi , multi , potential< multi >& , multi , int,multi, multi ,
     multi, multi, solution_grid< multi ,  std::vector<multi> , multi >&,
     multi&, multi, multi, multi, multi&,multi&, std::ostream&,bool,bool,
     multi,multi,bool,multi,multi,multi,bool,multi,multi,bool,bool);
//------------------------------------------------------------------------------
template DLL_EXPORT void odeSolve_dS_one_sided
    < multi , multi >
    (multi , multi , multi , potential< multi >& , multi , int,multi, multi ,
     multi, multi, solution_grid< multi ,  std::vector<multi> , multi >&,
     multi&, multi, multi, multi, multi&,multi&, std::ostream&,bool,bool);
//------------------------------------------------------------------------------
template DLL_EXPORT void odeSolve_dS_one_sided
    < multi , multi >
    (multi , multi , multi , potential< multi >& , multi , int,multi, multi ,
     multi, multi, solution_grid< multi ,  std::vector<multi> , multi >&,
     multi&, multi, multi, multi, multi&,multi&, std::ostream&,bool,bool,
     multi,multi,bool,multi);
//------------------------------------------------------------------------------
template DLL_EXPORT void odeSolve_dS_one_sided
    < multi , multi >
    (multi , multi , multi , potential< multi >& , multi , int,multi, multi ,
     multi, multi, solution_grid< multi ,  std::vector<multi> , multi >&,
     multi&, multi, multi, multi, multi&,multi&, std::ostream&,bool,bool,
     multi,multi,bool,multi,multi,multi,bool);
//------------------------------------------------------------------------------
/*
//------------------------------------------------------------------------------
//Double precision versions:
template DLL_EXPORT void odeSolve_dS
    < double , double >
    (double , double , double , potential< double >& , double , int,double,
     double ,
     double, double, solution_grid< double ,  std::vector<double> , double >&,
     solution_grid< double ,  std::vector<double> , double >&,
     solution_grid< double ,  std::vector<double> , double >&,
     double&, double, double, double,std::vector<multi>&,std::ostream&,bool,
     bool,
     double,double,bool,double,double,double,bool,double,double,double,double,
     bool,bool);
//------------------------------------------------------------------------------
template DLL_EXPORT void odeSolve_dSsingle
    < double , double >
    (double , double , double , double , potential< double >& , double , int,
      double,
     double , double, double, solution_grid< double , std::vector<double>
     , double >&,
     double&, double, double, double, std::ostream&,bool,bool);
//------------------------------------------------------------------------------
template DLL_EXPORT void odeSolve_dS_one_sided
    < double , double >
    (double , double , double , potential< double >& , double , int,double,
      multi ,
     double, double, solution_grid< double ,  std::vector<double> , double >&,
     double&, double, double, double, double&,double&, std::ostream&,bool,bool,
     double,double,bool,double,double,double,bool,double,double,bool);
//------------------------------------------------------------------------------
*/
