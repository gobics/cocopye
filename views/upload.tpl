<html>
    <head>
        <title>CoCoPyE Web</title>
    </head>
    <body>
        <h1 align="center">CoCoPyE Web</h1>

        <hr><br>

        <div align="center">
            <form action="/upload" method="post" enctype="multipart/form-data">
              <input type="file" name="fastaFile" /><br><br><br>
              <input type="submit" value="Submit" />
            </form>

            %if defined("result"):
                <br><h3>Results</h3>
                Completeness: {{ result[0] }}%<br>
                Contamination: {{ result[1] }}%
            %end

            %if defined("error"):
                <br><h3>Error</h3>
                {{ error }}
            %end
        </div>
    </body>
</html>