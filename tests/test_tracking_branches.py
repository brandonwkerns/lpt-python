import numpy as np
from context import lpt
import matplotlib.pylab as plt
from matplotlib import dates
import datetime as dt

plt.close('all')


options={}
options['min_overlap_points'] = 4000
options['min_overlap_frac'] = 0.5
options['min_lp_objects_points'] = 100
options['center_jump_max_hours'] = 3*24   # How long to allow center jumps
options['min_lpt_duration_hours'] = 7*24  # Minumum duration to keep it as an LPT

## Merging/Splitting settings
merge_split_options={}
merge_split_options['allow_merge_split'] = True
#merge_split_options['allow_merge_split'] = False
merge_split_options['split_merger_min_hours'] = 72     # Min duration of a split/merging track to separate it.

## For TMPA
options['objdir'] = '/home/orca/bkerns/public_html/realtime_mjo_tracking/lpt/data/tmpa/g20_72h/thresh12/objects'
#dt_list = [dt.datetime(2011,6,1,0,0,0) + dt.timedelta(hours=x) for x in np.arange(0.0, 396*24.0,3.0)]  #20*24
#dt_list = [dt.datetime(2011,11,1,0,0,0) + dt.timedelta(hours=x) for x in np.arange(0.0, 92*24.0,3.0)]  #20*24
dt_list = [dt.datetime(2011,12,8,0,0,0) + dt.timedelta(hours=x) for x in np.arange(0.0,20*24.0+1,3.0)]  #20*24
#dt_list = [dt.datetime(2011,12,20,0,0,0) + dt.timedelta(hours=x) for x in np.arange(0.0,9*24.0+1,3.0)]  #20*24


## For CMORPH
#options['objdir'] = '/home/orca/bkerns/public_html/realtime_mjo_tracking/lpt/data/cmorph/objects'
#dt_list = [dt.datetime(2019,3,1,0,0,0) + dt.timedelta(hours=x) for x in np.arange(0.0,31*24.0+1,1.0)]


print('Initialization.', flush=True)
LPT, BRANCHES = lpt.helpers.init_lpt_group_array(dt_list, options['objdir'])

print('Remove small LPOs.', flush=True)
LPT0, BRANCHES0 = lpt.helpers.lpt_group_array_remove_small_objects(LPT, BRANCHES, options)

print('Looping forward', flush=True)
LPTf, BRANCHESf = lpt.helpers.calc_lpt_group_array(LPT0, BRANCHES0, options, verbose=True)

print('Looping backwards', flush=True)
LPTb, BRANCHESb = lpt.helpers.calc_lpt_group_array(LPT0, BRANCHES0, options, verbose=True, reversed=True)

print('Merging forward and backwards', flush=True)
LPTfb0, BRANCHESfb0= lpt.helpers.overlap_forward_backward(LPTf, LPTb, BRANCHESf, BRANCHESb, options, verbose=True)

## Allow center jumps.
print(('Allow center jumps up to ' + str(options['center_jump_max_hours']) + ' hours.'), flush=True)
LPTfb1 = lpt.helpers.lpt_group_array_allow_center_jumps(LPTfb0, options)

print(('Remove LPT shorter than ' + str(options['min_lpt_duration_hours']) + ' hours.'), flush=True)
LPTfb2, BRANCHESfb2 = lpt.helpers.remove_short_lived_systems(LPTfb1, BRANCHESfb0, options['min_lpt_duration_hours'])

LPTfb, BRANCHESfb = lpt.helpers.lpt_group_id_separate_branches(LPTfb2, BRANCHESfb2, options, verbose=True)

print('Calc LP system group bulk properties.', flush=True)
TIMECLUSTERS0 = lpt.helpers.calc_lpt_system_group_properties_with_branches(LPTfb, BRANCHESfb, options)


fig0 = plt.figure(figsize=(50.0,30))
ax = fig0.add_subplot(121)
lpt.helpers.plot_timeclusters_time_lon(ax, TIMECLUSTERS0)
lpt.helpers.plot_lpt_groups_time_lon_text(ax, LPTfb, BRANCHESfb, options, text_color='k')
ax.set_title('Forward and Backwards, No Re-combinations')
ax.set_yticks(dt_list[::8])
ax.grid()
ax.yaxis.set_major_formatter(dates.DateFormatter("%m/%d"))


## Re-combine when splits/mergers are < 3 days.
if merge_split_options['allow_merge_split']:

    for this_group in np.unique(LPTfb[:,2]):
        print("!!!!!!!!!!!!!!!!!!!!!!!!!  Group #" + str(this_group) + "  !!!!!!!!!!!!!!!!!!!!!!!!!", flush=True)

        more_to_do = True
        niter=0

        """
        fig2 = plt.figure(figsize=(50.0,30))

        ax2 = fig2.add_subplot(111)
        lpt.helpers.plot_timeclusters_time_lon(ax2, TIMECLUSTERS0)
        lpt.helpers.plot_lpt_groups_time_lon_text(ax2, LPTfb, options, text_color='k')
        ax2.set_title('Iter # '+str(niter))
        ax2.set_yticks(dt_list[::8])
        ax2.grid()
        ax2.yaxis.set_major_formatter(dates.DateFormatter("%m/%d"))

        plt.savefig('test_tracking_iter'+str(niter).zfill(2)+'.png')

        plt.close(fig2)
        """




        ########################################################################
        ## Split and recombine cases. ##########################################
        ########################################################################
        print('------------------------------------------------')
        print('------------------------------------------------')
        print("Split and recombine (Rejoin 1).")
        print('------------------------------------------------')
        print('------------------------------------------------', flush=True)




        while more_to_do:
            more_to_do = False
            niter+=1

            #if niter > 5:
            #    break

            print('------------------------')
            print('Iteration #' + str(niter))
            print('------------------------', flush=True)

            branch_list = lpt.helpers.get_branches_in_lpt_group(LPTfb, BRANCHESfb, this_group)
            print("Unsorted branch list: " + str(branch_list))
            if len(branch_list) > 1:

                ## Put in order by duration.
                branch_durations = []

                for this_branch in branch_list:

                    lpt_all1 = lpt.helpers.lpt_branches_indices_list(LPTfb, BRANCHESfb, this_group, this_branch)
                    if len(lpt_all1) < 2:
                        continue
                    dt1all_list = sorted(np.unique([lpt.helpers.get_objid_datetime(LPTfb[xx,1]) for xx in lpt_all1]))
                    dt1all_begin = dt1all_list[0]
                    dt1all_end = dt1all_list[-1] #lpt.helpers.get_objid_datetime(LPTfb[np.max(lpt_diff1),1])

                    branch_durations += [(dt1all_end - dt1all_begin).total_seconds()/3600.0]

                branch_list_sorted = [branch_list[x] for x in np.argsort(branch_durations)]
                print("Sorted branch list: " + str(branch_list_sorted))

                for this_branch in branch_list_sorted:

                    #Which branch do you have the most intersection with?
                    other_branch = -1
                    max_intersect_duration = -1
                    for other_branch_test in branch_list_sorted:

                        if this_branch == other_branch_test:
                            continue

                        intersection_with_other_branch_test = lpt.helpers.lpt_branches_intersection(LPTfb,  BRANCHESfb, this_group, this_branch, other_branch_test)
                        if len(intersection_with_other_branch_test) > max_intersect_duration:
                            max_intersect_duration = len(intersection_with_other_branch_test)
                            other_branch = other_branch_test

                    if max_intersect_duration > 0:

                        print(str(this_branch) + ' with ' + str(other_branch) + '.')

                        lpt_all1 = lpt.helpers.lpt_branches_indices_list(LPTfb,  BRANCHESfb, this_group, this_branch)
                        if len(lpt_all1) < 2:
                            continue
                        dt1all_list = sorted(np.unique([lpt.helpers.get_objid_datetime(LPTfb[xx,1]) for xx in lpt_all1]))
                        dt1all_begin = dt1all_list[0]
                        dt1all_end = dt1all_list[-1] #lpt.helpers.get_objid_datetime(LPTfb[np.max(lpt_diff1),1])

                        lpt_all2 = lpt.helpers.lpt_branches_indices_list(LPTfb,  BRANCHESfb, this_group, other_branch)
                        if len(lpt_all2) < 2:
                            continue
                        dt2all_list = sorted(np.unique([lpt.helpers.get_objid_datetime(LPTfb[xx,1]) for xx in lpt_all2]))
                        dt2all_begin = dt2all_list[0]
                        dt2all_end = dt2all_list[-1] #lpt.helpers.get_objid_datetime(LPTfb[np.max(lpt_diff1),1])


                        lpt_diff1 = lpt.helpers.lpt_branches_difference(LPTfb, BRANCHESfb, this_group, this_branch, other_branch)
                        if len(lpt_diff1) > 0:

                            #print('Diff1: ' + str(lpt_diff1))
                            dt1_list = sorted(np.unique([lpt.helpers.get_objid_datetime(LPTfb[xx,1]) for xx in lpt_diff1]))
                            dt1_begin = dt1_list[0]
                            dt1_end = dt1_list[-1] #lpt.helpers.get_objid_datetime(LPTfb[np.max(lpt_diff1),1])
                            dur1 = (dt1_end - dt1_begin).total_seconds()/3600.0

                        else:

                            dt1_begin = None
                            dt1_end = None
                            dur1 = 0.0

                        lpt_diff2 = lpt.helpers.lpt_branches_difference(LPTfb, BRANCHESfb, this_group, other_branch, this_branch)
                        if len(lpt_diff2) > 0:
                            #print('Diff2: ' + str(lpt_diff2))
                            dt2_list = sorted(np.unique([lpt.helpers.get_objid_datetime(LPTfb[xx,1]) for xx in lpt_diff2]))
                            dt2_begin = dt2_list[0]
                            dt2_end = dt2_list[-1] #lpt.helpers.get_objid_datetime(LPTfb[np.max(lpt_diff1),1])
                            dur2 = (dt2_end - dt2_begin).total_seconds()/3600.0

                        else:

                            dt2_begin = None
                            dt2_end = None
                            dur2 = 0.0

                        ## Check to see if I have any difference in the branches.
                        ## NOTE: If branches have already been merged, this should not get triggered again.
                        if (len(lpt_diff1) > 0 or len(lpt_diff2) > 0):
                            ## If the difference is embeded in the intersection, it is a split-them-recombine case.
                            ## The two LPT branches are to be merged in to one.

                            if dt1_begin is None:
                                dt1_begin = dt2_begin #Note: datetimes are immutable.
                                dt1_end = dt2_end #Note: datetimes are immutable.
                            if dt2_begin is None:
                                dt2_begin = dt1_begin #Note: datetimes are immutable.
                                dt2_end = dt1_end #Note: datetimes are immutable.

                            #print((dt1_begin, dt1_end))
                            #print((dt2_begin, dt2_end))
                            #print((dt1all_begin, dt1all_end))
                            #print((dt2all_begin, dt2all_end))


                            ## Make sure the separated portion is *not* at the beginning or end of either track.
                            if (dt1_begin > dt1all_begin and dt1_end < dt1all_end) and (dt2_begin > dt2all_begin and dt2_end < dt2all_end):
                                print("Split and Re-combine.")
                                print("--> Combine these LPT branches.")


                                # Remove the smaller branch.
                                branches_to_remove = lpt.helpers.get_group_branches_as_list(BRANCHESfb[lpt_diff2[0]])
                                for branch_to_remove in branches_to_remove:
                                    LPTfb, BRANCHESfb = lpt.helpers.remove_branch_from_group(LPTfb, BRANCHESfb, this_group, branch_to_remove)

                                # Assign those LPT group array indices to the larger branch.
                                for jj in lpt_diff2:
                                    for kk in lpt_diff1:
                                        ## Only "inherit" larger branches for the relevant time span.
                                        kkdt = lpt.helpers.get_objid_datetime(LPTfb[kk,1])
                                        if kkdt >= dt2_begin and kkdt <= dt2_end:
                                            #LPTfb[jj,6] = int(LPTfb[jj,6]) | int(LPTfb[kk,6])
                                            BRANCHESfb[jj] = BRANCHESfb[jj] | BRANCHESfb[kk]

                                more_to_do = True

                                #print('Calc LP system group bulk properties.')
                                #TIMECLUSTERS0 = lpt.helpers.calc_lpt_system_group_properties_with_branches(LPTfb, options)

                                """
                                fig2 = plt.figure(figsize=(50.0,30))

                                ax2 = fig2.add_subplot(111)
                                lpt.helpers.plot_timeclusters_time_lon(ax2, TIMECLUSTERS0)
                                lpt.helpers.plot_lpt_groups_time_lon_text(ax2, LPTfb, options, text_color='k')
                                ax2.set_title('Iter # '+str(niter))
                                ax2.set_yticks(dt_list[::8])
                                ax2.grid()
                                ax2.yaxis.set_major_formatter(dates.DateFormatter("%m/%d"))

                                plt.savefig('test_tracking_iter'+str(niter).zfill(2)+'.png')

                                plt.close(fig2)

                                break
                                """

                    if more_to_do:
                        #more_to_do = False
                        break




        ########################################################################
        ## Splits and mergers. #################################################
        ########################################################################
        print('------------------------------------------------')
        print('------------------------------------------------')
        print('Splits and mergers (Rejoin 2).')
        print('------------------------------------------------')
        print('------------------------------------------------', flush=True)

        more_to_do = True

        while more_to_do:
            more_to_do = False
            niter+=1

            #if niter > 5:
            #    break

            print('------------------------')
            print('Iteration #' + str(niter))
            print('------------------------', flush=True)

            branch_list = lpt.helpers.get_branches_in_lpt_group(LPTfb, BRANCHESfb, this_group)
            print("Unsorted branch list: " + str(branch_list))
            if len(branch_list) > 1:

                ## Put in order by duration.
                branch_durations = []

                for this_branch in branch_list:

                    lpt_all1 = lpt.helpers.lpt_branches_indices_list(LPTfb,  BRANCHESfb, this_group, this_branch)
                    if len(lpt_all1) < 2:
                        continue
                    dt1all_list = sorted(np.unique([lpt.helpers.get_objid_datetime(LPTfb[xx,1]) for xx in lpt_all1]))
                    dt1all_begin = dt1all_list[0]
                    dt1all_end = dt1all_list[-1] #lpt.helpers.get_objid_datetime(LPTfb[np.max(lpt_diff1),1])

                    branch_durations += [(dt1all_end - dt1all_begin).total_seconds()/3600.0]

                branch_list_sorted = [branch_list[x] for x in np.argsort(branch_durations)]
                print("Sorted branch list: " + str(branch_list_sorted))

                for this_branch in branch_list_sorted:

                    #Which branch do you have the most intersection with?
                    other_branch = -1
                    max_intersect_duration = -1
                    for other_branch_test in branch_list_sorted:

                        if this_branch == other_branch_test:
                            continue

                        intersection_with_other_branch_test = lpt.helpers.lpt_branches_intersection(LPTfb,  BRANCHESfb, this_group, this_branch, other_branch_test)
                        if len(intersection_with_other_branch_test) > max_intersect_duration:
                            max_intersect_duration = len(intersection_with_other_branch_test)
                            other_branch = other_branch_test

                    if max_intersect_duration > 0:

                        print(str(this_branch) + ' with ' + str(other_branch) + '.')

                        lpt_all1 = lpt.helpers.lpt_branches_indices_list(LPTfb,  BRANCHESfb, this_group, this_branch)
                        if len(lpt_all1) < 2:
                            continue
                        dt1all_list = sorted(np.unique([lpt.helpers.get_objid_datetime(LPTfb[xx,1]) for xx in lpt_all1]))
                        dt1all_begin = dt1all_list[0]
                        dt1all_end = dt1all_list[-1] #lpt.helpers.get_objid_datetime(LPTfb[np.max(lpt_diff1),1])

                        lpt_all2 = lpt.helpers.lpt_branches_indices_list(LPTfb,  BRANCHESfb, this_group, other_branch)
                        if len(lpt_all2) < 2:
                            continue
                        dt2all_list = sorted(np.unique([lpt.helpers.get_objid_datetime(LPTfb[xx,1]) for xx in lpt_all2]))
                        dt2all_begin = dt2all_list[0]
                        dt2all_end = dt2all_list[-1] #lpt.helpers.get_objid_datetime(LPTfb[np.max(lpt_diff1),1])


                        lpt_diff1 = lpt.helpers.lpt_branches_difference(LPTfb, BRANCHESfb, this_group, this_branch, other_branch)
                        if len(lpt_diff1) > 0:

                            #print('Diff1: ' + str(lpt_diff1))
                            dt1_list = sorted(np.unique([lpt.helpers.get_objid_datetime(LPTfb[xx,1]) for xx in lpt_diff1]))
                            dt1_begin = dt1_list[0]
                            dt1_end = dt1_list[-1] #lpt.helpers.get_objid_datetime(LPTfb[np.max(lpt_diff1),1])
                            dur1 = (dt1_end - dt1_begin).total_seconds()/3600.0

                        else:

                            dt1_begin = None
                            dt1_end = None
                            dur1 = 0.0

                        lpt_diff2 = lpt.helpers.lpt_branches_difference(LPTfb, BRANCHESfb, this_group, other_branch, this_branch)
                        if len(lpt_diff2) > 0:
                            #print('Diff2: ' + str(lpt_diff2))
                            dt2_list = sorted(np.unique([lpt.helpers.get_objid_datetime(LPTfb[xx,1]) for xx in lpt_diff2]))
                            dt2_begin = dt2_list[0]
                            dt2_end = dt2_list[-1] #lpt.helpers.get_objid_datetime(LPTfb[np.max(lpt_diff1),1])
                            dur2 = (dt2_end - dt2_begin).total_seconds()/3600.0

                        else:

                            dt2_begin = None
                            dt2_end = None
                            dur2 = 0.0

                        ## Check to see if I have any difference in the branches.
                        ## NOTE: If branches have already been merged, this should not get triggered again.
                        if (len(lpt_diff1) > 0 or len(lpt_diff2) > 0):
                        #if (len(lpt_diff1) > 0 and len(lpt_diff2) > 0):

                            ## If the difference is embeded in the intersection, it is a split-them-recombine case.
                            ## The two LPT branches are to be merged in to one.

                            if dt1_begin is None:
                                dt1_begin = dt2_begin #Note: datetimes are immutable.
                                dt1_end = dt2_end #Note: datetimes are immutable.
                            if dt2_begin is None:
                                dt2_begin = dt1_begin #Note: datetimes are immutable.
                                dt2_end = dt1_end #Note: datetimes are immutable.

                            #print((dt1_begin, dt1_end))
                            #print((dt2_begin, dt2_end))
                            #print((dt1all_begin, dt1all_end))
                            #print((dt2all_begin, dt2all_end))


                            print("Merger or Split.")
                            #print(dur1)
                            #print(dur2)
                            if min(dur1,dur2) > merge_split_options['split_merger_min_hours'] - 0.1:
                                print("--> Retain both branches.")
                            else:
                                print("--> Combine these LPT branches.")
                                if dur1 > dur2:

                                    print('1 > 2')

                                    # Remove the smaller branch.
                                    branches_to_remove = lpt.helpers.get_group_branches_as_list(BRANCHESfb[lpt_diff2[0]])
                                    for branch_to_remove in branches_to_remove:
                                        LPTfb, BRANCHESfb = lpt.helpers.remove_branch_from_group(LPTfb, BRANCHESfb, this_group, branch_to_remove)

                                    # Assign those LPT group array indices to the larger branch.
                                    for jj in lpt_diff2:
                                        for kk in lpt_diff1:
                                            ## Only "inherit" larger branches for the relevant time span.
                                            kkdt = lpt.helpers.get_objid_datetime(LPTfb[kk,1])
                                            if kkdt >= dt2_begin and kkdt <= dt2_end:
                                                #LPTfb[jj,6] = int(LPTfb[jj,6]) | int(LPTfb[kk,6])
                                                BRANCHESfb[jj] = BRANCHESfb[jj] | BRANCHESfb[kk]

                                else:

                                    print('2 > 1')

                                    # Remove the smaller branch.
                                    branches_to_remove = lpt.helpers.get_group_branches_as_list(BRANCHESfb[lpt_diff1[0]])
                                    for branch_to_remove in branches_to_remove:
                                        LPTfb, BRANCHESfb = lpt.helpers.remove_branch_from_group(LPTfb, BRANCHESfb, this_group, branch_to_remove)

                                    # Assign those LPT group array indices to the larger branch.
                                    for jj in lpt_diff1:
                                        for kk in lpt_diff2:
                                            ## Only "inherit" larger branches for the relevant time span.
                                            kkdt = lpt.helpers.get_objid_datetime(LPTfb[kk,1])
                                            if kkdt >= dt1_begin and kkdt <= dt1_end:
                                                #LPTfb[jj,6] = int(LPTfb[jj,6]) | int(LPTfb[kk,6])
                                                BRANCHESfb[jj] = BRANCHESfb[jj] | BRANCHESfb[kk]

                                more_to_do = True


                                #print('Calc LP system group bulk properties.')
                                #TIMECLUSTERS0 = lpt.helpers.calc_lpt_system_group_properties_with_branches(LPTfb, options)

                                """
                                fig2 = plt.figure(figsize=(50.0,30))

                                ax2 = fig2.add_subplot(111)
                                lpt.helpers.plot_timeclusters_time_lon(ax2, TIMECLUSTERS0)
                                lpt.helpers.plot_lpt_groups_time_lon_text(ax2, LPTfb, options, text_color='k')
                                ax2.set_title('Iter # '+str(niter))
                                ax2.set_yticks(dt_list[::8])
                                ax2.grid()
                                ax2.yaxis.set_major_formatter(dates.DateFormatter("%m/%d"))

                                plt.savefig('test_tracking_iter'+str(niter).zfill(2)+'.png')

                                plt.close(fig2)
                                """
                                break


                    if more_to_do:
                        #more_to_do = False
                        break


print('Calc LP system group bulk properties.', flush=True)
TIMECLUSTERS0 = lpt.helpers.calc_lpt_system_group_properties_with_branches(LPTfb, BRANCHESfb, options)

ax = fig0.add_subplot(122)
lpt.helpers.plot_timeclusters_time_lon(ax, TIMECLUSTERS0)
lpt.helpers.plot_lpt_groups_time_lon_text(ax, LPTfb, BRANCHESfb, options, text_color='k')
ax.set_title('Forward and Backwards, With Recombinations')
ax.set_yticks(dt_list[::8])
ax.grid()
ax.yaxis.set_major_formatter(dates.DateFormatter("%m/%d"))

plt.savefig('test_tracking.png')

fn_tc = './TEST.txt'
lpt.lptio.lpt_system_tracks_output_ascii(fn_tc, TIMECLUSTERS0)

lpt.lptio.lpt_systems_group_array_output_ascii('TEST.group_array.txt', LPTfb, BRANCHESfb)


#plt.show()
