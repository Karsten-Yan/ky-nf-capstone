def getIndexes(dfObj, value):
    listOfPos = list()
    # Get bool dataframe with True at positions where the given value exists
    result = dfObj.isin([value])
    # Get list of columns that contains the value
    seriesObj = result.any()
    if seriesObj.any():
        columnNames = list(seriesObj[seriesObj == True].index)
        # Iterate over list of columns and fetch the rows indexes where value exists
        for col in columnNames:
            rows = list(result[col][result[col] == True].index)
        for row in rows:
            listOfPos.append((row))
        # Return a list of tuples indicating the positions of value in the dataframe
    else:
        listOfPos = False
    return listOfPos

def plot_dna(i, col, axs, ref, startpos, stoppos, df_temp):
    start,stop = 0,15
    font =  {'size': 24,
             'weight':'bold'}


    sequence = "".join(df_temp.base_3.values)

    modified_positions = getIndexes(df_temp[["base_3","modified_status"]], 1.0)

    if modified_positions:
            for elm in modified_positions:

                axs[0,col].axvspan(elm-0.5, elm+0.5, color='red', alpha=0.2,zorder=0)
                axs[1,col].axvspan(elm-0.5, elm+0.5, color='red', alpha=0.2,zorder=1)
                axs[2,col].axvspan(elm-0.5, elm+0.5, color='red', alpha=0.2,zorder=1)
                axs[3,col].axvspan(elm-0.5, elm+0.5, color='red', alpha=0.2,zorder=1)



    box1 = patches.Rectangle (xy=(4.5,1), width = 5, height = 2, edgecolor="grey",
                             linewidth=2.5, facecolor = "#FFB482")
    box2 = patches.FancyBboxPatch (xy=(6,5), width = 2, height = 2, edgecolor="grey",
                             linewidth=2.5, facecolor = "#A1C9F4", boxstyle = "Round4")

    axs[0,col].add_patch(box1)
    axs[0,col].add_patch(box2)
    axs[0,col].set_ylim(0,10)
    axs[0,col].text(7,2,s="reading frame", ha="center", va="center")
    axs[0,col].text(7,6,s="center", ha="center", va="center")
    axs[0,col].text(-0.3,2, va="center", ha="center", s="3'")
    axs[0,col].text(14.3,2, va="center", ha="center", s="5'")
    temp_title = df_temp.file_type[0]
    axs[0,col].text(7,9,s=temp_title, ha="center", va="center")
    axs[0,col].plot(np.arange(start,stop),np.full(stop,2), color="grey", zorder=0)
    axs[0,col].plot([4.5,6.5],[3,5],color = "grey", zorder=0)
    axs[0,col].plot([7.5,9.5],[5,3], color="grey", zorder=0)
    axs[0,col].axis("off")
    axs[0,col].set_yticklabels([""])

    bl = "#A1C9F4"
    gr = "#8DE5A1"
    og = "#FFB482"

    bar_colors = [gr,gr,gr,gr,gr,og,og,bl,og,og,gr,gr,gr,gr,gr,]

    axs[1,col].bar(x=range(0,stop-start),height=df_temp.dwell_time_diff_to_median, linewidth=2.5,
            edgecolor="grey",zorder=2, color = bar_colors)
    axs[1,col].set_ylabel("Dwell Time\nDiff to Median")
    axs[1,col].set_ylim(0,0.1)
    axs[1,col].set_yticklabels([""])

    axs[2,col].bar(x=range(0,stop-start),height=df_temp.dwell_time_rolling_min, linewidth=2.5,
            edgecolor="grey",zorder=2, color = bar_colors)
    axs[2,col].set_ylabel("Dwell Time\nRolling Minimum")
    axs[2,col].set_ylim(0,0.03)
    axs[2,col].set_yticklabels([""])

    axs[3,col].bar(x=range(0,stop-start),height=df_temp.dwell_time_median, linewidth=2.5,
            edgecolor="grey",zorder=2, color = bar_colors)
    axs[3,col].set_ylabel("Dwell Time\nMedian")
    axs[3,col].set_ylim(0,0.1)
    axs[3,col].set_yticklabels([""])
    axs[3,col].set_xticks(np.arange(start,stop))
    axs[3,col].set_xticklabels(list(sequence), fontdict=font )
    axs[3,col].set_xlim(-0.5,14.5)

    if df_temp.base_3.iloc[7] == "A":
        img=mpimg.imread('base_img/adenosin.png')
    elif df_temp.base_3.iloc[7] == "C":
        img=mpimg.imread('base_img/cytidin.png')
    elif df_temp.base_3.iloc[7] == "G":
        img=mpimg.imread('base_img/guanosin.png')
    elif df_temp.base_3.iloc[7] == "T":
        img=mpimg.imread('base_img/uridin.png')
    axs[1,2].imshow(img)
    axs[1,2].set_xlim(-100,1200)
    axs[1,2].set_ylim(1100,-100)
    axs[1,2].set_yticklabels([""])
    axs[1,2].set_xticklabels([""])
    axs[1,2].axis("off")

    if modified_positions:
            for elm in modified_positions:
                if elm == 7:
                    img2=mpimg.imread('base_img/1m7.png')
                    axs[2,2].imshow(img2)
    axs[2,2].set_xlim(-100,1400)
    axs[2,2].set_ylim(1100,-100)
    axs[2,2].set_yticklabels([""])
    axs[2,2].set_xticklabels([""])
    axs[2,2].axis("off")

def sequencing_show(df,i,ref):

    start_i = i
    stop_i =i+15

    df_temp_mod = df[(df["file_name"]=="modified_rep_1.tsv") &
                (df["ref_number"]=="ref_000"+ref)][start_i:stop_i].reset_index()
    df_temp_unmod = df[(df["file_name"]=="unmodified_rep_1.tsv") &
                (df["ref_number"]=="ref_000"+ref)][start_i:stop_i].reset_index()

    fig, axs = plt.subplots(
        4, 3,   figsize=(20, 10),
        sharex='col',
        gridspec_kw={"width_ratios": [2,2,1]}
    )


    plot_dna(i,0,axs,ref, start_i, stop_i, df_temp_mod)
    plot_dna(i,1,axs,ref, start_i, stop_i, df_temp_unmod)


    #plt.figure(edgecolor = "grey")
    plt.suptitle("DNA-Nanoporesequencing readout",weight="bold")
    #plt.tight_layout()
    plt.subplots_adjust(top=0.95,hspace = 0)
    axs[1,1].set_ylabel("")
    axs[2,1].set_ylabel("")
    axs[3,1].set_ylabel("")

    axs[0,2].set_visible(False)
    axs[3,2].bar(x = (400,900),
                 height= (df_temp_mod.dwell_time_median.iloc[7],df_temp_unmod.dwell_time_median.iloc[7]),
                 edgecolor="grey",
                 linewidth=2.5,
                 width=400)
    #axs[3,2].set_xlim(-0.5,1.5)
    axs[3,2].set_xticks([400,900])
    axs[3,2].set_xticklabels(["Mod","Unmod"] )
    axs[3,2].set_ylim(0,0.1)
    axs[3,2].set_yticklabels([""])
    axs[3,2].set_ylabel("Dwell Time\nMedian")
