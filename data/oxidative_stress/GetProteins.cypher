WITH "MATCH (p:protein {txid: 'txid224308'})-[r:ProGo]-(g:go_term {id: 'GO:0006979'}) RETURN p.id as id" AS query
CALL apoc.export.csv.query(query, "txid224308-stress-proteins.csv", {})
YIELD file, source, format, nodes, relationships, properties, time, rows, batchSize, batches, done, data
RETURN file, source, format, nodes, relationships, properties, time, rows, batchSize, batches, done, data;
WITH "MATCH (p:protein {txid: 'txid6239'})-[r:ProGo]-(g:go_term {id: 'GO:0006979'}) RETURN p.id as id" AS query
CALL apoc.export.csv.query(query, "txid6239-stress-proteins.csv", {})
YIELD file, source, format, nodes, relationships, properties, time, rows, batchSize, batches, done, data
RETURN file, source, format, nodes, relationships, properties, time, rows, batchSize, batches, done, data;
WITH "MATCH (p:protein {txid: 'txid7227'})-[r:ProGo]-(g:go_term {id: 'GO:0006979'}) RETURN p.id as id" AS query
CALL apoc.export.csv.query(query, "txid7227-stress-proteins.csv", {})
YIELD file, source, format, nodes, relationships, properties, time, rows, batchSize, batches, done, data
RETURN file, source, format, nodes, relationships, properties, time, rows, batchSize, batches, done, data;
WITH "MATCH (p:protein {txid: 'txid7955'})-[r:ProGo]-(g:go_term {id: 'GO:0006979'}) RETURN p.id as id" AS query
CALL apoc.export.csv.query(query, "txid7955-stress-proteins.csv", {})
YIELD file, source, format, nodes, relationships, properties, time, rows, batchSize, batches, done, data
RETURN file, source, format, nodes, relationships, properties, time, rows, batchSize, batches, done, data;
WITH "MATCH (p:protein {txid: 'txid559292'})-[r:ProGo]-(g:go_term {id: 'GO:0006979'}) RETURN p.id as id" AS query
CALL apoc.export.csv.query(query, "txid559292-stress-proteins.csv", {})
YIELD file, source, format, nodes, relationships, properties, time, rows, batchSize, batches, done, data
RETURN file, source, format, nodes, relationships, properties, time, rows, batchSize, batches, done, data;