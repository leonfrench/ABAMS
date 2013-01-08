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
package ubic.BAMSandAllen;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.AllenDataLoaders.AllenAtlasAnnotationLoader;
import ubic.BAMSandAllen.AllenDataLoaders.ImageSeriesInfoLoader;
import ubic.BAMSandAllen.AllenDataLoaders.StructureCatalogLoader;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes.Direction;
import ubic.BAMSandAllen.FocusedAnalysis.GeneCoordinateCorrelation;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenExpressionMatrixPair;
import ubic.BAMSandAllen.MatrixPairs.ExpressionMatrixPairFactory;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import ubic.basecode.dataStructure.params.ParamKeeper;
import ubic.basecode.math.CorrelationStats;
import ubic.basecode.math.DescriptiveWithMissing;
import ubic.basecode.math.Rank;
import ubic.basecode.math.Wilcoxon;
import cern.colt.list.DoubleArrayList;

public class TopTenInfo {
    private static Log log = LogFactory.getLog( TopTenInfo.class.getName() );

    RankedGeneListLoader loader;
    ConnectivityAndAllenExpressionMatrixPair pair;
    DoubleMatrix<String, String> expression;
    DoubleMatrix<String, String> connectivity;
    String baseFileName;
    Set<String> motorRegions;
    GeneCoordinateCorrelation forCoords;

    public TopTenInfo( RankedGeneListLoader loader, ConnectivityAndAllenExpressionMatrixPair pair ) throws Exception {
        super();
        this.baseFileName = loader.filename;
        this.loader = loader;
        this.pair = pair;
        // assume its a connectivity and expression pair, can be generalized later
        expression = pair.getMatrixB();
        connectivity = pair.getMatrixA();
        motorRegions = getMotorRegions();
        forCoords = new GeneCoordinateCorrelation( pair );

    }

    public void writeExpressionInfo() throws Exception {

//        AllenBrainAtlasService abaService = new AllenBrainAtlasService();
        ImageSeriesInfoLoader imageInfo = new ImageSeriesInfoLoader();

        getMotorRegions();

        List<String> rows;
        if ( loader != null ) {
            rows = loader.getRows();
        } else {
            // may be too many
            rows = expression.getRowNames();
        }

        DoubleArrayList averageExp = new DoubleArrayList();
        DoubleArrayList sdExp = new DoubleArrayList();
        for ( String geneRow : expression.getRowNames() ) {
            double[] expValues = expression.getRowByName( geneRow );

            DoubleArrayList expValuesDAL = new DoubleArrayList( expValues );
            double mean = DescriptiveWithMissing.mean( expValuesDAL );
            averageExp.add( mean );

            double sampleStandardDeviation = Math.sqrt( DescriptiveWithMissing.sampleVariance( expValuesDAL, mean ) );
            sdExp.add( sampleStandardDeviation );
        }
        DoubleArrayList averageRanks = Rank.rankTransform( averageExp );
        DoubleArrayList sdRanks = Rank.rankTransform( sdExp );

        log.info( "degree correlation to X:" + forCoords.getDegreeCorrelation( "x" ) );
        log.info( "degree correlation to Y:" + forCoords.getDegreeCorrelation( "y" ) );
        log.info( "degree correlation to Z:" + forCoords.getDegreeCorrelation( "z" ) );

        // //do motor regions have less connections
        DoubleArrayList motor1 = new DoubleArrayList();
        DoubleArrayList notMotor1 = new DoubleArrayList();

        ABAMSDataMatrix matrixA = pair.getMatrixA();
        DoubleMatrix<String, String> degNamed = Util.columnSums( matrixA );

        for ( String region : degNamed.getColNames() ) {
            double expvalue = degNamed.getByKeys( "Sums", region );
            if ( Double.isNaN( expvalue ) ) continue;
            if ( motorRegions.contains( region ) )
                motor1.add( expvalue );
            else
                notMotor1.add( expvalue );
        }
        log.info( "Wilcox Motor P " + Wilcoxon.exactWilcoxonP( motor1.elements(), notMotor1.elements() ) );
        log.info( "Wilcox Motor P Rev " + Wilcoxon.exactWilcoxonP( notMotor1.elements(), motor1.elements() ) );
        // //delete

        int rowsInMatrix = expression.rows();

        ParamKeeper stats = new ParamKeeper();
        int columns = pair.getMatrixA().columns();
        int count = 0;

        for ( String geneRow : rows ) {
            count++;
            if ( count % 100 == 0 ) log.info( "Count:" + count );

            Map<String, String> geneStats = new HashMap<String, String>();
            geneStats.put( "Name", geneRow );
            String geneName = ImageSeriesInfoLoader.getGeneNameFromRowName( geneRow );
            String imageSeriesID = ImageSeriesInfoLoader.getImageIDFromRowLabel( geneRow );
            // String imageSeriesID = geneRow.substring( geneRow.indexOf( "[" ) + 1, geneRow.indexOf( "]" ) );

            geneStats.put( "Gene", "'" + geneName );

            geneStats.put( "imageSeriesID", "'" + imageSeriesID );

            geneStats.put( "NCBI ID", "" + imageInfo.getNCBIIDFromRowName( geneRow ) );

            Set<String> allImageSets = imageInfo.getRowsFromGene( geneName );

            geneStats.put( "ImageSeriesCount", "" + allImageSets.size() );

            // do any of it's image series sets have a coronal image set?
            geneStats.put( "HasCoronalImage", "" + imageInfo.hasCoronalImageFromRowName( geneName ) );

            // find out how many other imageseries are in the top list
            allImageSets.retainAll( rows );
            geneStats.put( "ImageSeriesInList", "" + allImageSets.size() );

            // AbaGene gene = abaService.getGene( geneName );

            // geneStats.put( "plane", imageInfo.getPlaneFromRowName( geneRow ) );

            geneStats.put( "Gene Name", imageInfo.getNameFromImageID( imageSeriesID ) );

            // DoubleArrayList rx = Rank.rankTransform( x );

            int index = rows.indexOf( geneRow );
            geneStats.put( "Index", index + "" );

            int fullMatrixIndex = expression.getRowNames().indexOf( geneRow );
            geneStats.put( "meanRank", "" + ( rowsInMatrix - averageRanks.get( fullMatrixIndex ) ) );
            geneStats.put( "sdRank", "" + ( rowsInMatrix - sdRanks.get( fullMatrixIndex ) ) );

            double[] expValues = expression.getRowByName( geneRow );
            int nans = Util.countNaNs( expValues );
            geneStats.put( "NaNs", "" + nans );
            int columnsMinusNaN = columns - nans;

            double pearsonDegreeCor = pair.getDegreeCorrelation( geneRow );
            geneStats.put( "Degree correlation", pearsonDegreeCor + "" );
            geneStats.put( "Degree correlation pval", "" + CorrelationStats.pvalue( pearsonDegreeCor, columnsMinusNaN )
                    * rows.size() );

            double spearmanDegreeCor = pair.getRankDegreeCorrelation( geneRow );
            geneStats.put( "Rank degree correlation", spearmanDegreeCor + "" );
            geneStats.put( "Rank degree correlation pval", ""
                    + CorrelationStats.spearmanPvalue( spearmanDegreeCor, columnsMinusNaN ) * rows.size() );

            boolean removeNaNs = true;
            geneStats.put( "expSum", "" + Util.sum( expValues, removeNaNs ) );

            geneStats.put( "plane", imageInfo.getPlaneFromRowName( geneRow ) );

            double coordCorrelation = forCoords.getXRankCorrelation( geneRow );
            geneStats.put( "x correlation", "" + coordCorrelation );
            geneStats.put( "x correlation pvalue", ""
                    + CorrelationStats.spearmanPvalue( coordCorrelation, columnsMinusNaN ) * rows.size() );

            coordCorrelation = forCoords.getYRankCorrelation( geneRow );
            geneStats.put( "y correlation", "" + coordCorrelation );
            geneStats.put( "y correlation pvalue", ""
                    + CorrelationStats.spearmanPvalue( coordCorrelation, columnsMinusNaN ) * rows.size() );
            coordCorrelation = forCoords.getZRankCorrelation( geneRow );
            geneStats.put( "z correlation", "" + coordCorrelation );
            geneStats.put( "z correlation pvalue", ""
                    + CorrelationStats.spearmanPvalue( coordCorrelation, columnsMinusNaN ) * rows.size() );

            DoubleArrayList expValuesDAL = new DoubleArrayList( expValues );
            double mean = DescriptiveWithMissing.mean( expValuesDAL );
            geneStats.put( "mean", "" + mean );

            double sampleStandardDeviation = Math.sqrt( DescriptiveWithMissing.sampleVariance( expValuesDAL, mean ) );
            geneStats.put( "sampleStandardDeviation", "" + sampleStandardDeviation );

            // // do test for genes enriched in motor regions
            DoubleArrayList motor = new DoubleArrayList();
            DoubleArrayList notMotor = new DoubleArrayList();
            for ( String region : expression.getColNames() ) {
                double expvalue = expression.getByKeys( geneRow, region );
                if ( Double.isNaN( expvalue ) ) continue;
                if ( motorRegions.contains( region ) )
                    motor.add( expvalue );
                else
                    notMotor.add( expvalue );
            }
            geneStats.put( "Wilcox Motor P", "" + Wilcoxon.exactWilcoxonP( motor.elements(), notMotor.elements() )
                    * rows.size() );
            geneStats.put( "Wilcox Motor P Rev", "" + Wilcoxon.exactWilcoxonP( notMotor.elements(), motor.elements() )
                    * rows.size() );

            // url
            String url = "http://mouse.brain-map.org/brain/" + geneName + "/" + imageSeriesID + ".html?ispopup=true";
            url = "HYPERLINK(\"" + url + "\",\"Series summary\")";
            geneStats.put( "URL", url );

            String urlThumbs = "http://mouse.brain-map.org/brain/" + geneName + "/" + imageSeriesID
                    + "/thumbnails.html?ispopup=true";
            urlThumbs = "HYPERLINK(\"" + urlThumbs + "\",\"Thumbnails\")";
            geneStats.put( "Thumbs URL", urlThumbs );

            // regions above average?

            stats.addParamInstance( geneStats );

        }
        stats.writeExcel( baseFileName + ".geneInfo.xls" );
    }

    private Set<String> getMotorRegions() throws Exception {
        StructureCatalogLoader allenCatalog = new StructureCatalogLoader();
        Set<String> motorRegionsAllen = allenCatalog.getMotorNuclei();
        Set<String> motorRegions = new HashSet<String>();

        // convert to BAMS names
        for ( String allenRegion : motorRegionsAllen ) {
            motorRegions.addAll( pair.convertBNametoA( allenRegion ) );
        }
        return motorRegions;
    }

    public void writeColumnInfo() throws Exception {
        ParamKeeper stats = new ParamKeeper();

        for ( String colname : connectivity.getColNames() ) {
            Map<String, String> colStats = getColInfo( colname );
            stats.addParamInstance( colStats );
        }

        stats.writeExcel( baseFileName + ".colInfo.xls" );
    }

    private Map<String, String> getColInfo( String colname ) throws Exception {
        Map<String, String> colStats = new HashMap<String, String>();

        colStats.put( "Name", colname );
        double[] conValues = connectivity.getColumnByName( colname );
        String connections = "";
        for ( String conRow : connectivity.getRowNames() ) {
            if ( connectivity.getByKeys( conRow, colname ) == 1d ) {
                connections += "|" + conRow;
            }
        }
        colStats.put( "connections", connections );

        double[] expValues = expression.getColumnByName( colname );
        colStats.put( "connectionCount", "" + Util.sum( conValues ) );

        // assume one allen name
        StructureCatalogLoader allenCatalog = new StructureCatalogLoader();

        colStats.put( "MotorRegion", "" + motorRegions.contains( colname ) );
        colStats.put( "x", "" + forCoords.getX( colname ) );
        colStats.put( "y", "" + forCoords.getY( colname ) );
        colStats.put( "z", "" + forCoords.getZ( colname ) );

        if ( pair.isVirtualRegion( colname ) ) {
            colStats.put( "allen names", colname );
            colStats.put( "virtual", "true" );
        } else {
            colStats.put( "virtual", "false" );
            Set<String> allennames = allenCatalog.getAllenMappedRegions( colname );
            String allenName;
            colStats.put( "allen names", allennames.toString() );
            if ( !allennames.isEmpty() ) {
                allenName = allennames.iterator().next();
                AllenAtlasAnnotationLoader spaceLoader = new AllenAtlasAnnotationLoader();
                DoubleMatrix<String, String> volumeMatrix = spaceLoader.getVolumeMatrix();
                colStats.put( "volume", "" + volumeMatrix.getByKeys( "volume", allenName ) );
            }
        }

        boolean removeNaNs = true;
        colStats.put( "expSum", "" + Util.sum( expValues, removeNaNs ) );
        colStats.put( "expNaNs", "" + Util.countNaNs( expValues ) );
        return colStats;
    }

    /**
     * @param args
     */
    public static void main( String[] args ) throws Exception {
        boolean removeNonExp = true;
        boolean useVirtual = false;
        Direction direction = AnalyzeBAMSandAllenGenes.Direction.APPENDED;
        ConnectivityAndAllenExpressionMatrixPair pair = ExpressionMatrixPairFactory.connectivityAndExpression(
                direction, useVirtual, removeNonExp );

//        log.info( pair.getCorrelation() );
//        pair.writeRMatrices();
//        System.exit( 1 );
        //
        // pair.shuffleDataRows( 1 );

        RankedGeneListLoader loader = new RankedGeneListLoader( pair.getMatrixBDataRows(), "AllGenesAppendedInfo" );

        loader = new RankedGeneListLoader(
                "/grp/java/workspace/BAMSandAllen/data/ABACellEnriched/patterns/top100-pattern-b.txt", false );

        TopTenInfo test = new TopTenInfo( loader, pair );
        // test.writeColumnInfo();
        test.writeExpressionInfo();
        System.exit( 1 );
        //
        // // String endFile = "LOOGenesInOrder.outgoing.partialcon.txt";
        // String endFile = "LOOGenesInOrder.outgoing.partialcon.txt.topGenes.txt";
        // String baseFolder = SetupParameters.config.getString( "abams.dataFolder" ) + "rankedGenes/near final/";
        // // String outputFile = SetupParameters.config.getString( "abams.dataFolder" ) + "/topten/" + endFile;
        // String outputFile = baseFolder + endFile;

        // RankedGeneListLoader loader = new RankedGeneListLoader( baseFolder + endFile, false );

        // TopTenInfo test = new TopTenInfo( loader, pair );
        // test.writeColumnInfo();
        // // test.writeExpressionInfo();
    }
}
