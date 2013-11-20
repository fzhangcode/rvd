//
//  main.c
//  pileup2dc
//
//  Created by Patrick Flaherty on 3/27/12.
//  Copyright (c) 2012 Stanford University. All rights reserved.
//

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>

#define MAXLINE  5000000 // maximum line length
#define MAXDEPTH 2500000 // maximum depth
enum base { A,C,G,T,N };

// data structure for a single position
typedef struct {
    int pos; // 1-base position in reference
    enum base refb; // reference base
    
    int depthF[5], depthR[5], depthT[5]; // forward/reverse depth in base enum order
    
} pile_t;

// get pileup file length
static int get_filelength ( FILE *fid ) {
    int numLines = 0;
    char *tempString = (char *) malloc (MAXLINE * sizeof(char));
    // need to use malloc because tempString[BIGNUM] causes stack overflow
 
    rewind(fid); // start counting from the file start
    while( !feof(fid) ) {
        fgets(tempString, MAXLINE, fid);
        numLines++;
    }
    rewind(fid);
    

    free (tempString);
    return numLines-1;
}

// compute total depth
static int sum_depth(const int *d) {
    return d[0]+d[1]+d[2]+d[3]+d[4];
}

// parse a single alignment string
static void parse_align_string(const char refb, const char *alignString, pile_t *pile) 
{
    char alignChar;
    int i, indelLen;
    
    
    switch (refb) {
        case 'A':
            pile->refb=A;
            break;
        case 'C':
            pile->refb=C;
            break;
        case 'G':
            pile->refb=G;
            break;
        case 'T':
            pile->refb=T;
            break;   
        case 'N':
            pile->refb=N;
            break; 
        default:
            break;
    }
    
    int alignLen = (int)strlen(alignString);
    
    
    for (i=0; i<alignLen; i++) {
        alignChar = alignString[i];
        switch (alignChar) {
            case '.': // forward reference
                pile->depthF[pile->refb]++;
                break;
            case ',': //reverse reference
                pile->depthR[pile->refb]++;
                break;
                
                
            case 'A': 
                pile->depthF[A]++;
                break;
            case 'C':
                pile->depthF[C]++;
                break;
            case 'G':
                pile->depthF[G]++;
                break;
            case 'T':
                pile->depthF[T]++;
                break;
            case 'N':
                pile->depthF[N]++;
                break;
                
                // reverse read different from reference
            case 'a':
                pile->depthR[A]++;
                break;
            case 'c':
                pile->depthR[C]++;
                break;
            case 'g':
                pile->depthR[G]++;
                break;
            case 't':
                pile->depthR[T]++;
                break;
            case 'n':
                pile->depthR[N]++;
                break;
                
            // TODO: Is it a good idea to skip indels?
            case '+':
            case '-': //skip indels
                i++;
                indelLen = alignString[i]-'0';
                i = i+indelLen;
                break;
            
            case '^': // ignore start of read segment and the following mapping quality
                i++;
                break;
            case '$': // ignore the end of the read segment
                break;
                
            default: // ignore * characters. TODO: where are these coming from?
                break;
        }
    }
    
    // accumulate the forward and reverse depths
    for(i = 0; i<5; i++) {
        pile->depthT[i] = pile->depthF[i]+pile->depthR[i];
    }
    
    return;
}

// get the majority vote base
static int get_majvot (const int *depth)
{
    int i;
    
    int majVot = 0, majDepth=0;
    for (i = 0; i<5; i++) {
        if (depth[i] > majDepth) {
            majVot = i;
            majDepth = depth[i];
        }
    }
    return majVot;
}


// export the data in depth chart format
static void export_depthchart(const pile_t *pile, const int numLines, const char *chr)
{
    //char *chr="chr1";
    char baseChar[5]="ACGTN";
    int i;
    
    // Export header information
    fprintf(stdout, "count\tchr\tloc\tsubref\trefb");
    fprintf(stdout, "\t");
    
    fprintf(stdout, "A_All\tC_All\tG_All\tT_All\tN_All\tTot_All\tmatch_All\tMajVot_All");
    fprintf(stdout, "\t");
    
    fprintf(stdout, "A_F\tC_F\tG_F\tT_F\tN_F\tTot_F\tmatch_F\tMajVot_F");
    fprintf(stdout, "\t");
    
    fprintf(stdout, "A_R\tC_R\tG_R\tT_R\tN_R\tTot_R\tmatch_R\tMajVot_R");
    fprintf(stdout, "\n");
    
    for (i=0; i<numLines; i++) {
        fprintf(stdout, "%d\t%s\t%d\t%s\t%c", i+1, chr, pile[i].pos, chr, baseChar[pile[i].refb]);
        fprintf(stdout, "\t");
        fprintf(stdout, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%c", pile[i].depthT[A], pile[i].depthT[C], 
                pile[i].depthT[G], pile[i].depthT[T], pile[i].depthT[N],
                sum_depth(pile[i].depthT), pile[i].depthT[pile[i].refb], baseChar[get_majvot(pile[i].depthT)]);
        fprintf(stdout, "\t");
        fprintf(stdout, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%c", pile[i].depthF[A], pile[i].depthF[C], 
                pile[i].depthF[G], pile[i].depthF[T], pile[i].depthF[N],
                sum_depth(pile[i].depthF), pile[i].depthF[pile[i].refb], baseChar[get_majvot(pile[i].depthF)]);
        fprintf(stdout, "\t");
        fprintf(stdout, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%c", pile[i].depthR[A], pile[i].depthR[C], 
                pile[i].depthR[G], pile[i].depthR[T], pile[i].depthR[N],
                sum_depth(pile[i].depthR), pile[i].depthR[pile[i].refb], baseChar[get_majvot(pile[i].depthR)]);
        fprintf(stdout, "\n");
    }
}

static void init_depth(pile_t *pile)
{
    int i;

    for(i=0; i<5; i++) { 
        pile->depthF[i] = 0;
        pile->depthR[i] = 0;
        pile->depthT[i] = 0;
    }
    
}

int main (int argc, char *argv[])
{
    char *lineBuf = (char *) malloc (MAXLINE * sizeof(char));
    char chrName[40], refb;
    char *alignString = (char *) malloc (MAXDEPTH * sizeof(char));
    char *qualityString = (char *) malloc (MAXDEPTH * sizeof(char));
    int pos, totDepth;
    int i;
    FILE *plfid;
    
    if (argc == 1) {
        fprintf(stderr, "Usage: pileup2dc [in.pileup]\n");
        return 1;
    }
    
    
    // open the pileup file for reading
    // die if error opening or reading file
    plfid = fopen(argv[1], "rd");
 
    // get the number of positions in the pileup
    int numLines = get_filelength(plfid);

    // allocate space for the depth chart
    pile_t pile[numLines];
    
    // read the pileup file and store in the structure
    for (i=0; i<numLines; i++) {

		fgets(lineBuf, MAXLINE, plfid);
        sscanf(lineBuf, "%s\t%d\t%c\t%d", chrName, &pos, &refb, &totDepth);

		init_depth(&pile[i]);
		pile[i].pos = pos;
    

		
		if ( totDepth > 0) {
			sscanf(lineBuf, "%s\t%d\t%c\t%d\t%s\t%s", chrName, &pos, &refb,
               &totDepth, alignString, qualityString);
			parse_align_string ( toupper(refb), alignString, &pile[i] );
		}
		else {
			switch (refb) {
        			case 'A':
            			pile[i].refb=A;
            			break;
        			case 'C':
            			pile[i].refb=C;
            			break;
        			case 'G':
            			pile[i].refb=G;
            			break;
        			case 'T':
            			pile[i].refb=T;
            			break;   
        			case 'N':
            			pile[i].refb=N;
           			break; 
        			default:
            			break;
    			}
		}
    }
    fclose(plfid);
    
    // export the pileup structure in depth chart format
    export_depthchart ( pile, numLines, chrName);
    
    free(lineBuf);
    free(alignString);
    free(qualityString);
    return 0;
}
