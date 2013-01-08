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

import java.util.HashSet;
import java.util.Set;
import java.util.concurrent.CopyOnWriteArraySet;

import org.apache.commons.lang.time.StopWatch;

import ubic.BAMSandAllen.ClassSelectors.ClassSelector;
import ubic.basecode.dataStructure.CountingMap;

import com.hp.hpl.jena.ontology.OntClass;
import com.hp.hpl.jena.ontology.OntModel;
import com.hp.hpl.jena.rdf.model.Property;
import com.hp.hpl.jena.rdf.model.RDFNode;
import com.hp.hpl.jena.rdf.model.Resource;
import com.hp.hpl.jena.rdf.model.Statement;
import com.hp.hpl.jena.rdf.model.StmtIterator;
import com.hp.hpl.jena.util.iterator.ExtendedIterator;
import com.hp.hpl.jena.vocabulary.RDFS;

public class JenaUtil {

    public static Set<OntClass> fitlerClasses( ExtendedIterator ei, ClassSelector filter ) {
        return fitlerClasses( ei.toSet(), filter );
    }

    public static Set<OntClass> fitlerClasses( OntModel m, ClassSelector filter ) {
        return fitlerClasses( m.listClasses(), filter );
    }

    public static Set<OntClass> fitlerClasses( Set classes, ClassSelector filter ) {
        Set<OntClass> result = new CopyOnWriteArraySet<OntClass>();
        for ( Object o : classes ) {
            OntClass ontClass = ( OntClass ) o;
            if ( filter.test( ontClass ) ) {
                result.add( ontClass );
            }
        }
        return result;
    }

    public static Set<Resource> getSubjects( StmtIterator it ) {
        Set<Resource> result = new HashSet<Resource>();
        while ( it.hasNext() ) {
            result.add( it.nextStatement().getSubject() );
        }
        return result;
    }

    public static CountingMap<String> getLiteralStringCounts( StmtIterator it ) {
        CountingMap<String> literals = new CountingMap<String>();
        while ( it.hasNext() ) {
            literals.increment( it.nextStatement().getLiteral().getLexicalForm() );
        }
        return literals;
    }

    public static Set<Resource> getObjects( StmtIterator it ) {
        Set<Resource> result = new HashSet<Resource>();
        while ( it.hasNext() ) {
            RDFNode object = it.nextStatement().getObject();
            if ( object.isResource() ) {
                result.add( ( Resource ) ( object.as( Resource.class ) ) );
            }
        }
        return result;
    }

    public static String getLabel( Resource r ) {
        return getStringLiteral( r, RDFS.label );
    }

    public static String getStringLiteral( Resource r, Property p ) {
        Statement s = r.getProperty( p );
        if ( s == null ) return null;
        return s.getLiteral().getLexicalForm();
    }
}
