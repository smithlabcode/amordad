-- Created by Vertabelo (http://vertabelo.com)
-- Last modification date: 2015-04-04 17:57:06.997




-- tables
-- Table feature_vector
CREATE TABLE feature_vector (
    id varchar(255)    NOT NULL ,
    path varchar(255)    NOT NULL ,
    CONSTRAINT feature_vector_pk PRIMARY KEY (id)
) ENGINE=INNODB;

-- Table graph_edge
CREATE TABLE graph_edge (
    src varchar(255)    NOT NULL ,
    dst varchar(255)    NOT NULL ,
    dist double(10,9)    NOT NULL ,
    CONSTRAINT graph_edge_pk PRIMARY KEY (src,dst)
) ENGINE=INNODB;

-- Table hash_function
CREATE TABLE hash_function (
    id varchar(255)    NOT NULL ,
    path varchar(255)    NOT NULL ,
    update_time timestamp    NOT NULL ,
    CONSTRAINT hash_function_pk PRIMARY KEY (id)
) ENGINE=INNODB;

-- Table hash_table_bucket
CREATE TABLE hash_table_bucket (
    id varchar(255)    NOT NULL ,
    hash_key int    NOT NULL ,
    occupant varchar(255)    NOT NULL ,
    CONSTRAINT hash_table_bucket_pk PRIMARY KEY (id,hash_key,occupant)
) ENGINE=INNODB;





-- foreign keys
-- Reference:  destination_node (table: graph_edge)


ALTER TABLE graph_edge ADD CONSTRAINT destination_node FOREIGN KEY destination_node (dst)
    REFERENCES feature_vector (id)
    ON DELETE CASCADE
    ON UPDATE CASCADE;
-- Reference:  hash_table (table: hash_table_bucket)


ALTER TABLE hash_table_bucket ADD CONSTRAINT hash_table FOREIGN KEY hash_table (id)
    REFERENCES hash_function (id)
    ON DELETE CASCADE
    ON UPDATE CASCADE;
-- Reference:  occupant (table: hash_table_bucket)


ALTER TABLE hash_table_bucket ADD CONSTRAINT occupant FOREIGN KEY occupant (occupant)
    REFERENCES feature_vector (id)
    ON DELETE CASCADE
    ON UPDATE CASCADE;
-- Reference:  source_node (table: graph_edge)


ALTER TABLE graph_edge ADD CONSTRAINT source_node FOREIGN KEY source_node (src)
    REFERENCES feature_vector (id)
    ON DELETE CASCADE
    ON UPDATE CASCADE;



-- End of file.

