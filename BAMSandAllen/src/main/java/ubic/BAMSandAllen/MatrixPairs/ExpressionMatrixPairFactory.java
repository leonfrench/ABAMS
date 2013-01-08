/*******************************************************************************
 * The ABAMS project
 * 
 * Copyright (c) 2012 University of British Columbia
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *       http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 ******************************************************************************/
package ubic.BAMSandAllen.MatrixPairs;

import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.ABAMSDataMatrix;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes.Direction;
import ubic.BAMSandAllen.ClassSelectors.BrainRegionClassSelector;
import ubic.BAMSandAllen.FocusedAnalysis.ExploreRegionNames;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenPartialExpressionMatrixPair.RegressMatrix;
import ubic.BAMSandAllen.geneFilters.ABA1000ConflictsGeneFilter;
import ubic.BAMSandAllen.geneFilters.AGEAGeneFilter;
import ubic.BAMSandAllen.geneFilters.GeneFilter;
import ubic.BAMSandAllen.geneFilters.NaNGeneFilter;
import ubic.BAMSandAllen.geneFilters.NonExpFilter;
import ubic.BAMSandAllen.geneFilters.NullNoneGeneFilter;
import ubic.basecode.dataStructure.StringToStringSetMap;

public class ExpressionMatrixPairFactory {
    private static Log log = LogFactory.getLog( ExpressionMatrixPairFactory.class.getName() );

    public static ABAMSDataMatrix getEnergyMatrix( Direction direction, boolean removeNonExp ) throws Exception {
        // don't run, leave brain regions as Allen
        boolean run = false;
        boolean useVirtual = false;
        boolean square = false;
        ConnectivityAndAllenExpressionMatrixPair pair = connectivityAndExpression( direction, useVirtual, removeNonExp,
                run, square );
        return pair.getMatrixB();
    }

    public static ConnectivityAndAllenExpressionMatrixPair connectivityAndExpression( Direction direction,
            boolean useVirtual, boolean removeNonExp ) throws Exception {
        boolean run = true;
        boolean square = false;
        return connectivityAndExpression( direction, useVirtual, removeNonExp, run, square );
    }

    public static ConnectivityAndAllenExpressionMatrixPair connectivityAndExpression( Direction direction,
            boolean useVirtual, boolean removeNonExp, boolean run, boolean square ) throws Exception {
        long t;
        double zeroReplace = Double.NaN;
        boolean log1p = false;
        boolean doLog = true;
        ConnectivityAndAllenExpressionMatrixPair pair = new ConnectivityAndAllenExpressionMatrixPair(
                new BrainRegionClassSelector(), doLog, log1p, square, zeroReplace, "NewEnergies", direction );

        if ( useVirtual ) {
            pair.makeVirtualRegions();
        }

        if ( removeNonExp ) {
            pair.applyGeneFilter( new NonExpFilter() );
        }

        if ( run ) {
            pair.run();
        }

        // default filters
        for ( GeneFilter filter : ExpressionMatrixPairFactory.getDefaultExpressionFilters() ) {
            pair.applyGeneFilter( filter );
        }

        if ( run ) {
            log.info( pair.getCorrelation() );
        }

        return pair;
    }

    public static ConnectivityAndAllenPartialExpressionMatrixPair connectivityPartial( Direction direction,
            boolean slow, RegressMatrix regressType, boolean useVirtual, boolean removeNonExp, boolean logDistance )
            throws Exception {
        boolean square = false;
        return connectivityPartial( direction, slow, regressType, useVirtual, removeNonExp, logDistance, square );
    }

    public static ConnectivityAndAllenPartialExpressionMatrixPair connectivityPartial( Direction direction,
            boolean slow, RegressMatrix regressType, boolean useVirtual, boolean removeNonExp, boolean logDistance,
            boolean square ) throws Exception {
        boolean run = true;
        return connectivityPartial( direction, slow, regressType, useVirtual, removeNonExp, logDistance, square, run );
    }

    public static ConnectivityAndAllenPartialExpressionMatrixPair connectivityPartial( Direction direction,
            boolean slow, RegressMatrix regressType, boolean useVirtual, boolean removeNonExp, boolean logDistance,
            boolean square, boolean run ) throws Exception {
        if ( regressType == RegressMatrix.CONNECTIVITY ) {
            slow = false;
            log.info( "Setting to fast mode for regression on connectivity" );
        }
        // this pair is used to convert the space pair region names to BAMS names
        ConnectivityAndAllenSpacePair spacePair = new ConnectivityAndAllenSpacePair( new BrainRegionClassSelector(),
                square, null, direction, logDistance );
        if ( useVirtual ) spacePair.makeVirtualRegions();
        spacePair.run();

        ConnectivityAndAllenPartialExpressionMatrixPair partialMantelPair = new ConnectivityAndAllenPartialExpressionMatrixPair(
                new BrainRegionClassSelector(), true, false, square, Double.NaN, "NewEnergies", direction, spacePair
                        .getMatrixB(), regressType );

        // make virtual regions
        if ( useVirtual ) partialMantelPair.makeVirtualRegions();

        if ( removeNonExp ) {
            partialMantelPair.applyGeneFilter( new NonExpFilter() );
        }

        // remove a region??

        if ( run ) partialMantelPair.run();

        // default filters
        for ( GeneFilter filter : ExpressionMatrixPairFactory.getDefaultExpressionFilters() ) {
            partialMantelPair.applyGeneFilter( filter );
        }

        if ( run ) log.info( partialMantelPair.getCorrelation() );
        return partialMantelPair;
    }

    public static List<String> getUsedCols( boolean squareMatrix, Direction direction ) throws Exception {
        // most of these settings (log, log1, log replace) don't matter as its only getting the columns
        // don't run same space as that will convert it to BAMS names
        ConnectivityAndAllenExpressionMatrixPair expressionPair = new ConnectivityAndAllenExpressionMatrixPair(
                new BrainRegionClassSelector(), true, false, squareMatrix, Double.NaN, "NewEnergies", direction );
        // forR.writeRMatrices();
        expressionPair.slimMatrices();
        expressionPair.printDimensions();
        expressionPair.removeConnectionZeroes();
        return expressionPair.getAllenDataColNames();
    }

    public static List<GeneFilter> getDefaultExpressionFilters() {
        return filters;
    }

    static List<GeneFilter> filters = new LinkedList<GeneFilter>();

    public static void addFilter( GeneFilter f ) {
        filters.add( f );
    }

    static {
        filters.add( new NullNoneGeneFilter() );
        filters.add( new NaNGeneFilter() );
        // NaN filter?
        // low var filter
    }
}
