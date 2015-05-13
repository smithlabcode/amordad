-- sql statement to make sure each node has no more than max_deg neighbors
select feature_vector.id, count(graph_edge.dst) as OutDegree
       from graph_edge left join feature_vector
       on feature_vector.id = graph_edge.src group by src;
