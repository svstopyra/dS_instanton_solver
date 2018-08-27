//Backend for the instanton finder in the AdS/Flat spacetime case.
#define DLL_EXPORT
#define INSTANTON_SOLVER_DLL_API //Comment this out if we wish to import these
//from a DLL instead.

#include "project_specific.h"

#include "multi_precision_definitions.h"
#include <boost/array.hpp>//For boost arrays (needed for ode_types)
#include "AdS.h"
#include <initializer_list>//Allow us to use initialiser lists for vector.
    //This simplifies the code somewhat.


//------------------------------------------------------------------------------
//Main function to export:
template< class value_type , class time_type >
void odeSolveAdSFlat(value_type false_vacuum , value_type true_vacuum ,
                     value_type barrier , potential< value_type >& V ,
                     time_type chimax , int odeSolverToUse,
                     value_type RelTol, value_type AbsTol, value_type stepError,
                     value_type xi,
                     value_type lowerBound , value_type upperBound ,
                     solution_grid< time_type , std::vector< value_type > ,
                     value_type >& solOut,value_type& DSout,
                     value_type precision,std::ostream& outStream,
                     bool track_scale_factor,bool compute_H0_linearisation)
{
//INPUTS:
    //Ok. we proceed by bisecting until we have a solution which overshoots
    //beyond
    //chimax, which stands in for infinity.

    //Overshoots occur when the solution crosses zero before chimax.
    //Undershoots when
    //its derivative crosses zero before then.

    //Instantiate the events function, state-type, dynamic BCs etc...
    typedef std::vector< value_type > AdSFlatSolType;//solution without action,
        //for speed.
    typedef std::vector< value_type > AdSFlatSolTypeAction;
    //typedef odeAdSFlat< AdSFlatSolType , time_type , value_type>odeRHSAdSFlat;
    typedef odeAdSFlat_delta_a< AdSFlatSolType , time_type , value_type>
            odeRHSAdSFlat;
    typedef eventsAdSFlat< AdSFlatSolType , time_type , value_type > evAdSFlat;
    /*typedef dynamicBCsAdSFlat_taylor< value_type , AdSFlatSolType , time_type>
            dynBCAdSFlat;*/
    typedef dynamicBCs_flatfv_delta< value_type , AdSFlatSolType , time_type>
            dynBCAdSFlat;
    /*typedef odeAdSFlatAction<AdSFlatSolTypeAction,time_type,value_type>
            odeRHSAdSFlatAction;*/
    typedef odeAdSFlatAction_delta_a<AdSFlatSolTypeAction,time_type,value_type>
            odeRHSAdSFlatAction;
    typedef eventsAdSFlat< AdSFlatSolTypeAction , time_type , value_type >
            evAdSFlatAction;
    typedef dynamicBCsAdSFlat_taylor< value_type , AdSFlatSolTypeAction ,
                                      time_type > dynBCAdSFlatAction;
    typedef eventsGravFriction<AdSFlatSolType,time_type,value_type> eventsGrav;
    typedef odeGravFriction< AdSFlatSolType , time_type , value_type > odeGrav;
    typedef dynamicBCsGravFriction_taylor< value_type , AdSFlatSolType ,
                                      time_type > dynBCGrav;
    typedef odeGravFrictionAction<AdSFlatSolType,time_type,value_type>
                                    odeGravAction;

    //Instantiation for overshoot/undershoot:

    const value_type _1p0 = value_type(1.0);
    const value_type _0p0 = value_type(0.0);
    const value_type _2p0 = value_type(2.0);
    const value_type _6p0 = value_type(6.0);
    std::cout.precision(50);

    //Potential
    value_type W0 = value_type(0.0);
    value_type h;
    if(xi >= _0p0 && xi < _1p0/_6p0)
    {
        h = _1p0;//h = phi_scale/M_p so h = 1 is Planck units.
    }
    else
    {
        h = _1p0/sqrt(abs(xi));
        //h = value_type(1.0);//h = phi_scale/M_p so h = 1 is Planck units.
    }
	value_type V0 = _0p0;//V at phi = 0.
    odeRHSAdSFlat odeToSolveTrack(V,xi, h, V0); //ode to solve (gravitational
                                                           //instanton equation)
    odeGrav odeToSolveNoTrack(V,xi, h, V0);
    odeRHS<std::vector<value_type>,time_type,value_type>* odeToSolve;
    if(track_scale_factor)
    {
        odeToSolve = &odeToSolveTrack;
    }
    else
    {
        odeToSolve = &odeToSolveNoTrack;
    }

    evAdSFlat eventsTrack(false_vacuum,true_vacuum,false_vacuum,outStream);
    //events function.
    eventsGrav eventsNoTrack(false_vacuum,true_vacuum);
    events_function< std::vector<value_type> , time_type, value_type >* events;
    if(track_scale_factor)
    {
        events = &eventsTrack;
    }
    else
    {
        events = &eventsNoTrack;
    }
        //Checks for events labelled IE:
        //IE = 1 -> y crosses zero (terminal event, => overshoot)
        //IE = 2 -> y' crosses zero (terminal event, => undershoot)
        //IE = 3 -> a crosses zero (non-terminal event, for reference only)
        //IE = 4 -> y crosses lowerBound (usually 0, but could be different -
                                            //terminal event)
        //IE = 5 -> y crosses upperBound (terminal event, indicates bad initial
                                          //condition or perhaps a->0)
    dynBCAdSFlat bcAdjusterTrack(xi,h,V0,V,stepError,1);
    dynBCGrav bcAdjusterNoTrack(xi,h,V0,V,stepError,1);
    //Performs analytic
    // first
    //step using a taylor series approximation.
        //stepError - determines size of initial step.
        //version (last argument) - specifies which version of the ode we are
            //using.
            //version = 0 -> normal ode, no action.
            //version = 1 -> ode with action
            //version = 2 -> linearised ode.
    dynamicBCs< value_type , std::vector<value_type> , time_type >* bcAdjuster;
    if(track_scale_factor)
    {
        bcAdjuster = &bcAdjusterTrack;
    }
    else
    {
        bcAdjuster = &bcAdjusterNoTrack;
    }






    //Search loop:
    int nMax = 200;
    int counter = 0;
    //Current bounds on location of bounce:
    value_type lower = lowerBound;//Barrier
    value_type upper = upperBound;//True vacuum
    outStream << "Lower initial = " << lower << std::endl;
    outStream << "Upper initial = " << upper << std::endl;
    value_type guess;
    multi diff = abs(upper - lower);
    outStream << "diff = " << diff << std::endl;
    while(diff > precision)
    {
        //outStream << "Loop stage 1" <<std::endl;
        //Bisect to get initial guess:
        guess = (lower + upper)/_2p0;
        outStream << "\nn = " << counter << ", guess = " << guess << "\n";
        outStream << "\ndiff = " << upper - lower;
        //outStream << "Loop stage 2" <<std::endl;

        //Search for FWHMif(track_scale_factor)
        if(track_scale_factor)
        {
            eventsTrack.value_to_check = (guess + false_vacuum)/_2p0;
        }



        //Formulate initial conditions vector:
        AdSFlatSolType y0;
        //value_type y0List[] = { guess , _0p0 , _0p0 , _1p0 };
        value_type y0List[] = { guess , _0p0 , _0p0 , _1p0 , _0p0, _0p0 };

        y0.assign(y0List,y0List + 6);
        //AdSFlatSolType y0NoTrack;
        //value_type y0NoTrackList[] = { guess , _0p0 , _0p0 , _1p0};
        //y0NoTrack.assign(y0TrackList,y0TrackList + 4);

        //AdSFlatSolType y0 = track_scale_factor ? y0Track : y0NoTrack;
        //Step away to avoid co-ordinate singularity at chi = 0:
        //outStream << "Loop stage 3" <<std::endl;
        AdSFlatSolType y01step = y0;
        time_type epsilon = (*bcAdjuster)(y0,y01step);//Copies solution after
            //step of epsilon into y01step. Last argument tells us whether
            //we step forwards (0) or backwards (1) in time.
        //Set up time range:
        //outStream << "Loop stage 4" <<std::endl;
        //time_type tspanArray[2] = { epsilon , chimax };
        std::vector< time_type > tspan;
        time_type tspan_list[] = { epsilon , chimax };
        tspan.assign(tspan_list,tspan_list + 2);
        //tspan.assign(tspanArray,tspanArray + 2);
        //Solve the ode:
        //outStream << "Loop stage 5" <<std::endl;
        solution_grid< time_type , AdSFlatSolType , value_type > solTemp;
        time_type initStepMax = time_type(0.01)/time_type(chimax);
        //outStream << "Integrating..." <<std::endl;
        int nSuccess = ode_solve<time_type,AdSFlatSolType,value_type>
            (*odeToSolve,y01step,tspan,(*events),odeSolverToUse,solTemp,RelTol,
             AbsTol,initStepMax);
        if(nSuccess != 1)
        {
            outStream << "nSuccess = " << nSuccess
                      << " when solving for AdS-Flat instantons." << std::endl;
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
                if(track_scale_factor)
                {
                    overshoots = true;
                    upper = guess;
                }
                //else we just found a' = 0, so not terminal.
                break;
            case 2:
                //undershoot, IF it occurs on the opposite side of the barrier
                //to the true vacuum (theoretically impossible to occur on the
                //other side, so must be a numerical artefact if it does).
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
                //Found a = 0 point.
                outStream << "Found a = 0. Error, or possibly background is "
                          << "de-Sitter like." << std::endl;
                //reportToCaller("Found a = 0. Error, or possibly background is
                                 //de-Sitter like.");
                throw "Unexpected ode-event encountered during integration.";
            case 4:
                //crossed upper bound, so overshoot for tracking a'/a version,
                //FWHM for the a',a separately tracked version.
                if(!track_scale_factor)
                {
                    overshoots = true;
                    upper = guess;
                    eventsNoTrack.y0upperBound = upper;
                }
                break;
            case 5:
                //Crossed lower bound (in both cases).
                //Usually indicates an overshoot:
                //outStream << "Overshoot." << std::endl;
                //reportToCaller("Overshoot.\n");
                overshoots = true;
                upper = guess;
                if(!track_scale_factor)
                {
                    eventsNoTrack.y0upperBound = upper;
                }
                break;
            case 6:
                //Crossed upper bound. Only
                //occurs for the a'/a version.
                //Indicates a bad initial condition.
                outStream << "Upper bound crossed during ode integration. "
                          << "Possibly a bad initial condition (wrong side of "
                          << "true vacuum) or solution undershot and event "
                          << "detection failed to detect this." << std::endl;
                //reportToCaller("Upper bound crossed during ode integration.
                    //Possibly a bad initial condition (wrong side
                    //of true vacuum) or solution undershot and event detection
                    //failed to detect this.\n");
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
        //outStream << "y0Upper - y0lower = " << diff << std::endl;
        if(overshoots)
        {
            outStream << "Overshoot.\n" << std::endl;
        }
        else if(undershoots)
        {
            outStream << "Undershoot.\n" << std::endl;
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
                for(int j = 0;j < y01step.size();j++)
                {
                    outStream << solTemp.YE[i][j] << std::endl;
                }
                outStream << "TE[" << i << "] = " << solTemp.TE[i] << std::endl;
            }
            outStream << "y0 = \n";
            for(int i = 0;i < y01step.size();i++)
            {
                outStream << y01step[i] << std::endl;
            }
            outStream << "yend = \n";
            for(int i = 0;i < y01step.size();i++)
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
    odeRHSAdSFlatAction odeToSolveActionTrack(V,xi,h,V0);
    evAdSFlatAction eventsActionTrack(false_vacuum,true_vacuum,false_vacuum,
                                 outStream);
    dynBCAdSFlatAction bcAdjusterActionTrack(xi,h,V0,V,stepError,2);
    odeGravAction odeToSolveActionNoTrack(V,xi,h,V0);
    dynBCGrav bcAdjusterActionNoTrack(xi,h,V0,V,stepError,2);
    odeRHS<std::vector<value_type>,time_type,value_type>* odeToSolveAction;
    dynamicBCs< value_type , std::vector<value_type> , time_type >*
        bcAdjusterAction;
    eventsGrav eventsActionNoTrack(false_vacuum,true_vacuum);
    events_function< std::vector<value_type> , time_type, value_type >*
        eventsAction;
    if(track_scale_factor)
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
    //outStream << "Computing final bounce solution... " << std::endl;

    //Check for when the solution passes through its FWHM:
    //This feature isn't yet implemented for the NoTrack version, so we don't
    //use it in that case.
    if(track_scale_factor)
    {
        eventsActionTrack.value_to_check = (guess + false_vacuum)/_2p0;
    }

    //Boundary conditions at x = 0
    AdSFlatSolTypeAction y0Action;
    //value_type y0ActionList[] = {guess,_0p0,_0p0,_1p0,_0p0};
    value_type y0ActionList[] = {guess,_0p0,_0p0,_1p0,_0p0,_0p0,_0p0};
    //y0Action.assign(y0ActionList,y0ActionList + 5);
    y0Action.assign(y0ActionList,y0ActionList + 7);
    AdSFlatSolTypeAction y0Action1step = y0Action;
    //Same BCs work for NoTrack, where we use theta = log(a + 1)

    //Step to boundary conditions slightly offset from x = 0:
    time_type epsAction = (*bcAdjusterAction)(y0Action,y0Action1step);

    //Integration range:
    time_type tspanArray[2] = { epsAction , chimax };
    std::vector< time_type > tspan;
    tspan.assign(tspanArray,tspanArray + 2);
    time_type initStepMax = time_type(0.01)/time_type(chimax);

    //Integrate:
    int nSuccess = ode_solve<time_type,AdSFlatSolTypeAction,value_type>
    (*odeToSolveAction,y0Action1step,tspan,*eventsAction,odeSolverToUse,solOut,
     RelTol,AbsTol,initStepMax);
    if(nSuccess != 1)
    {
        outStream << "nSuccess = " << nSuccess << " when solving for AdS-Flat"
                  << "instantons." << std::endl;
        throw "Integration failure.";
    }
    DSout = solOut.Y[int(solOut.T.size()) - 1][4];
    outStream << "Done." << std::endl;
    //Done.
}
//------------------------------------------------------------------------------
template< class value_type , class time_type >
void odeSolveAdSFlatSingle(value_type y0, value_type false_vacuum ,
                           value_type true_vacuum , value_type barrier ,
                           potential< value_type >& V , time_type chimax ,
                           int odeSolverToUse,
                           value_type RelTol, value_type AbsTol,
                           value_type stepError, value_type xi,
                           value_type lowerBound , value_type upperBound ,
                           solution_grid< time_type , std::vector< value_type >,
                           value_type >& solOut,value_type& DSout,
                           value_type precision,std::ostream& outStream,
                           bool track_scale_factor)
{
//INPUTS:

     outStream << "Starting Integration for y0 = \n" << y0 << std::endl;
     outStream << "\nV(y0) = " << V(y0) << std::endl;

     const value_type _1p0 = value_type(1.0);
     const value_type _0p0 = value_type(0.0);
     const value_type _2p0 = value_type(2.0);
     const value_type _6p0 = value_type(6.0);

    //Ok. we proceed by bisecting until we have a solution which overshoots
    //beyond
    //chimax, which stands in for infinity.

    //Overshoots occur when the solution crosses zero before chimax.
    //Undershoots when
    //its derivative crosses zero before then.

    //Instantiate the events function, state-type, dynamic BCs etc...
    typedef std::vector< value_type > AdSFlatSolType;//solution without action,
        //for speed.
    typedef std::vector< value_type > AdSFlatSolTypeAction;
    typedef odeAdSFlat< AdSFlatSolType , time_type , value_type > odeRHSAdSFlat;
    typedef eventsAdSFlat< AdSFlatSolType , time_type , value_type > evAdSFlat;
    typedef dynamicBCsAdSFlat_taylor< value_type , AdSFlatSolType , time_type >
            dynBCAdSFlat;
    typedef odeAdSFlatAction< AdSFlatSolTypeAction , time_type , value_type >
            odeRHSAdSFlatAction;
    typedef eventsAdSFlat< AdSFlatSolTypeAction , time_type , value_type >
            evAdSFlatAction;
    typedef dynamicBCsAdSFlat_taylor< value_type , AdSFlatSolTypeAction ,
                                      time_type > dynBCAdSFlatAction;
    typedef eventsGravFriction<AdSFlatSolType,time_type,value_type> eventsGrav;
    typedef odeGravFriction< AdSFlatSolType , time_type , value_type > odeGrav;
    typedef dynamicBCsGravFriction_taylor< value_type , AdSFlatSolType ,
                                      time_type > dynBCGrav;
    typedef odeGravFrictionAction<AdSFlatSolType,time_type,value_type>
                                    odeGravAction;

    //Instantiation for overshoot/undershoot:

    //Compute action:
    value_type W0 = _0p0;
    value_type h;
    if(xi >= _0p0 && xi < _1p0/_6p0)
    {
        h = _1p0;//h = phi_scale/M_p so h = 1 is Planck units.
    }
    else
    {
        h = _1p0/sqrt(abs(xi));
        //h = value_type(1.0);//h = phi_scale/M_p so h = 1 is Planck units.
    }
    //Setup odeRHS, event-detector, and dynamicBC objects needed:
    odeRHSAdSFlatAction odeToSolveActionTrack(V,xi,h,W0);
    evAdSFlatAction eventsActionTrack(false_vacuum,true_vacuum,false_vacuum,
                                 outStream);
    dynBCAdSFlatAction bcAdjusterActionTrack(xi,h,W0,V,stepError,2);
    odeGravAction odeToSolveActionNoTrack(V,xi,h,W0);
    dynBCGrav bcAdjusterActionNoTrack(xi,h,W0,V,stepError,2);
    odeRHS<std::vector<value_type>,time_type,value_type>* odeToSolveAction;
    dynamicBCs< value_type , std::vector<value_type> , time_type >*
        bcAdjusterAction;
    eventsGrav eventsActionNoTrack(false_vacuum,true_vacuum);
    events_function< std::vector<value_type> , time_type, value_type >*
        eventsAction;
    if(track_scale_factor)
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


    AdSFlatSolTypeAction y0Action;
    value_type y0ActionList[] = { y0 , _0p0 , _0p0 , _1p0 , _0p0 };
    y0Action.assign(y0ActionList,y0ActionList + 5);
    AdSFlatSolTypeAction y0Action1step = y0Action;
    outStream << "Compute boundary conditions.\n";
    time_type epsAction = (*bcAdjusterAction)(y0Action,y0Action1step);
    outStream << "Done.";
    outStream << "\ny0 = " << y0Action1step[0];
    outStream << "\ny0p = " << y0Action1step[1];
    outStream << "\na0 = " << y0Action1step[2];
    outStream << "\na0p = " << y0Action1step[3];
    outStream << "\nS = " << y0Action1step[4];
    time_type tspanArray[2] = { epsAction , chimax };
    std::vector< time_type > tspan;
    tspan.assign(tspanArray,tspanArray + 2);
    time_type initStepMax = time_type(0.01)/time_type(chimax);
    outStream << "Starting integration.";
    int nSuccess = ode_solve<time_type,AdSFlatSolTypeAction,value_type>
    (*odeToSolveAction,y0Action1step,tspan,*eventsAction,odeSolverToUse,solOut,
     RelTol,AbsTol,initStepMax);
    if(nSuccess != 1)
    {
        outStream << "nSuccess = " << nSuccess << " when solving for AdS-Flat"
                  << " instantons." << std::endl;
        throw "Integration failure.";
    }
    DSout = solOut.Y[int(solOut.T.size()) - 1][4];
    //Done.
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//Explicit instantiation and export of this function:
template DLL_EXPORT void odeSolveAdSFlat
    < multi , multi >
    (multi , multi , multi , potential< multi >& , multi , int,multi, multi ,
     multi, multi, multi , multi ,
     solution_grid< multi , std::vector< multi > , multi >&,
     multi&, multi, std::ostream&,bool,bool);
//------------------------------------------------------------------------------
template DLL_EXPORT void odeSolveAdSFlatSingle
    < multi , multi >
    (multi , multi , multi , multi , potential< multi >& , multi , int, multi,
     multi , multi, multi, multi , multi ,
     solution_grid< multi , std::vector< multi > , multi >&,
     multi&, multi, std::ostream&,bool);
//------------------------------------------------------------------------------
