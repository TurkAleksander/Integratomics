#This code contains chunks of the HelloWorld script for integratomics
#I am putting it here so it doesn't bother me as I work on the local version

class HelloWorld:
	
	def ping(self, **kw):
		cherrypy.session.load()
		sid = cherrypy.session.id
		
		if cherrypy.session.has_key('studyweights'):
			perm = int(cherrypy.session.get('studyweights')['perm'])
			interval = int(cherrypy.session.get('studyweights')['interval'])
			print ("Number of permutations is %s" % perm)
		else:
			cherrypy.session['perm'] = "1"
			cherrypy.session['name'] = "Ales"
			cherrypy.session.save()
		
		render = web.template.render('/home/ales/integr/website/templates/')
		yield render.testComet()['__body__']

		if cherrypy.session.has_key('perm'):
			a = int(cherrypy.session.get('perm')) + 10
			cherrypy.session['perm'] = str(a)
		else:
			cherrypy.session['perm'] = "1"
		yield '''<img style="margin-left:auto;margin-right:auto;display:block;" src="07p-6606-loadingcircle.gif" />'''
		yield '''<div id="sliders2">'''

		#yield '<script type="text/javascript">window.scrollBy(0,50);</script>'
		#yield "Session is %s" % cherrypy.session.get('perm')
		#yield cherrypy.session.get('rez1')['Maver']
		#integration(100, 100000, [1,1])
		interval=int(interval)
		perm = int(perm)
		overlap = interval/2
		np.set_printoptions(threshold='nan')
		
		#################################################################
		######## IMPORT DATA FROM ANNOTATION DATABASES
		#################################################################

		importAnnot()

		#################################################################
		######## Prepare genomic backbone based on chromosomal lengths
		#################################################################
				
		numRegions = prepareBackbone(interval, overlap)
		idx = numRegions

		#################################################################
		######## IMPORT DATA FOR GENE DENSITY
		#################################################################
		
		genes = hgncImport(interval, overlap, numRegions)		
		#################################################################
		######## IMPORT DATA FROM STUDIES
		#################################################################

		studyFiles = []
		studyTypes = []
		studyNames = []
		annoTypes=[]

		# Get weights for study types
		studyloc = '/home/ales/integr/temp/%s/Study*.txt' % sid
		for name in glob.glob(studyloc):
			z=open(name, 'r')
			studyTypes.append(re.sub(r'\W+', '', z.readline().rstrip().replace('-', '').replace(' ', '')))
			studyNames.append(re.sub(r'\W+', '', z.readline().rstrip().replace('-', '').replace(' ', '')))
			annoTypes.append(re.sub(r'\W+', '', z.readline().rstrip()))
			z.close()	
			
		# Detect syudyTypes
		uniqStudyTypes = unique(studyTypes)
		typeW = cherrypy.session.get('typeW')

		# List files in the directory
		for name in glob.glob(studyloc):
			studyFiles.append(name)
		
		numStudies = len(studyFiles)	
		mat=np.zeros((idx,numStudies))
		
		con = lite.connect('/home/ales/integr/sqlite3/test.db')
		cur = con.cursor()
		
		#################################################################
		######## PUT THE DATA INTO A MATRIX
		#################################################################		
		
		for study in range(0,numStudies):
			if annoTypes[study] == "coor":
				#Coordinate information reading
				output = importCoorStudy(studyFiles, study, studyTypes, studyNames, interval, overlap, cur, con)
								
			else:
				#Annotated information reading
				output = importAnnotStudy(studyFiles, study, studyTypes, studyNames, annoTypes, interval, overlap, cur, con)

			yield "Cleaning the data... <br />"
			output = clearZeros(np.reshape(output, idx))
			yield "Data cleaned! <br />"


			yield """##############################<br />"""
			yield "%s successfully imported into the database!<br />" % studyNames[study]
			yield """##############################<br />"""

			mat[:,study]=output

		#################################################################
		######## DO PERMUTATIONAL TEST
		#################################################################

		studyW = np.array([1]*numStudies)
		print (studyW)

		studyTypes = np.array(studyTypes)
		mat2 = np.zeros((numRegions, len(unique(studyTypes))))

		
		for i in range(len(unique(studyTypes))):
			print (i)
			type=unique(studyTypes)[i]
			perMat1 = perMat(mat[:,studyTypes==type], studyW[studyTypes==type], perm, genes)
			
			#perMat1.getSimpleScores()
			#mat2[:,i]=perMat1.simplescores
			
			for permNo in range(perm):
				perMat1.permutationTest()
				yield "Permutation %s... DONE!" % permNo
				yield '<br />'
			perMat1.getScores()
			mat2[:,i]=np.array(perMat1.scores)

		studyW = [1]*len(unique(studyTypes))
		perMat2 = perMat(mat2, typeW, perm, genes)

		for permNo in range(perm):
			perMat2.permutationTest()
			yield "Permutation %s... DONE!" % permNo
			yield '<br />'
		perMat2.getScores()
		
		np.savetxt("/home/ales/integr/perm_scores.txt", -np.log10(perMat2.scores))
		
		#################################################################
		######## SORT RESULTS
		#################################################################

		ntop=5000
		resultsIdx = np.arange(perMat2.nrow)[np.argsort(perMat2.scores)[::-1]][np.arange(ntop)]
		resultsScores = np.array(perMat2.scores[np.argsort(perMat2.scores)[::-1]][np.arange(ntop)])

		#Save a matrix with data for each study	
		np.savetxt("/home/ales/integr/myfile.txt", mat)

		#convert getscores to make an array with scores
		con = lite.connect('/home/ales/integr/sqlite3/test.db')
		cur = con.cursor()

		sql0 = "DROP TABLE IF EXISTS %s" % "Results"
		cur.execute(sql0)
		sql1 = "CREATE TABLE %s (idx INT, score INT)" % "Results"
		cur.execute(sql1)
		cur.execute("ALTER TABLE Backbone ADD COLUMN Results INT")
		sql2 = "UPDATE Backbone SET Results=? WHERE idx=?"
		print (resultsIdx.tolist())
		cur.executemany(sql2, zip(resultsScores.tolist(),resultsIdx.tolist()))
		cur.execute("SELECT Genes.name, max(Backbone.chr), max(Backbone.stop), min(Backbone.start), max(Backbone.Results) FROM Genes INNER JOIN Backbone ON Genes.start = Backbone.start AND Genes.chr=Backbone.chr AND Results>2 GROUP BY Genes.name")
		#sql3 = "SELECT * FROM Gene LEFT OUTER JOIN %s ON %s.start = Backbone.start AND Backbone.chr = %s.chr GROUP BY Backbone.idx" % (studyNames[study], studyNames[study], studyNames[study], studyNames[study])
		#cur.execute(sql3)
		con.commit()
		output = cur.fetchall()
		#print output

		try:
			tdir = '/home/ales/integr/temp/%s' % sid
			os.mkdir(tdir)
		except:
			print ("Folder already there!")

		rfile = '/home/ales/integr/temp/%s/results.txt' % sid
		f = open(rfile, 'w')
		for out in output:
			print (out)
			f.write(out[0] + '\t' + out[1] + '\t' + str(out[3]) + '\t' + str(out[2]) + '\t' + str(round(out[4],1)) + '\n') 
		f.close()
		
		webresultFile = '/home/ales/integr/website/static/uploadify/wgplots/%s.results.txt' % sid
		shutil.copyfile(rfile, webresultFile)
		
		cur.execute("SELECT idx FROM Backbone")
		output = cur.fetchall()
		idxs = np.reshape(np.array(output, dtype='int'), numRegions)
		#print idxs
		
		cur.execute("SELECT chr FROM Backbone")
		output = cur.fetchall()
		chrs = np.reshape(np.array(output, dtype='str'), numRegions)

		
		GWplot1 = GWplot(idxs, chrs, perMat2.scores, sid)
		
		yield '''<br />'''
		yield '''################# COMPLETED ##################<br />'''
		yield '''<br />'''
		yield 'You will be redirected to results page in '
		time.sleep(1)
		yield '3......   '
		time.sleep(1)
		yield '2......   ' 
		time.sleep(1)
		yield '1......   '
		yield '</div></body></html>'
		cherrypy.session.save()
		yield '''<script>window.location = '/results';</script>'''
		
	ping.exposed = True
	ping._cp_config = {'response.stream': True}

	def results(self):
		cherrypy.session.load()
		sid = cherrypy.session.id
		render = web.template.render('/home/ales/integr/website/templates/')
		yield render.results("Ales")['__body__']
		yield '''<div style="max-width:100%;"><img src="wgplots/''' + sid + '''.png"></div>'''
		
		webresultFile = 'wgplots/%s.results.txt' % sid
		yield '''<div style="max-width:100%;"><a href="''' + webresultFile + '''">Download results directly HERE</a></div>'''

		yield '''</div> 
			<div id="sliders1">
				<form>
					Filter by score:     <label class="cutoff-label-sliding" for="cutoff-sliding">Cutoff</label>
					<input type="text" name="cutoff" class="cutoff-sliding" id="cutoff" />
				</form><br />
			<span id="foo"></span>'''

	results.exposed = True
	
	def generateresults(self, cutoff):
		cherrypy.session.load()
		sid = cherrypy.session.id
		try:
			print ("Cutoff exists %s" % cutoff)
			float(cutoff)
		except:
			cutoff = 0
			print ("Cutoff doesnt exist %s" % cutoff)	
		try:
			float(cutoff)
		except:	
			return 'Please enter a numeric value in the cutoff field...'
						
		rfile = '/home/ales/integr/temp/%s/results.txt' % sid
		f = open(rfile, 'r')
		out=""
		for line in f.readlines():
			cols = line.rstrip().split('\t')
			if float(cols[4]) > float(cutoff):
				print ("OK!")
				out += '<tr><td>'+'</td><td>'.join([str(i) for i in cols])+'</td></tr>'
			else: 
				next
		return '''<style>#center{margin-left: auto;margin-right: auto;}</style>'''+'''<p>Top resulting genes are presented below</p>'''+'''<table id="table-3">
			<thead>
				<th>Gene Name</th>
				<th>Chromosome Name</th>
				<th>Region start</th>
				<th>Region stop</th>
				<th>Score</th>
			</thead>
			<tbody>'''+out+'</table>'
	generateresults.exposed = True		
	
	def index1(self):
		if cherrypy.session.has_key('perm'):
			a = int(cherrypy.session.get('perm')) + 1
			cherrypy.session['perm'] = str(a)
		else:
			cherrypy.session['perm'] = "1"
		cherrypy.session['name'] = "Ales"
		cherrypy.session.save()
		render = web.template.render('/home/ales/integr/website/templates/')
		yield render.uploadify("Ales")['__body__']
	index1.exposed = True
	
	def studyweights(self, **kw):
		cherrypy.session.load()
		sid = cherrypy.session.id		
		if cherrypy.session.has_key('perm'):
			a = int(cherrypy.session.get('perm')) + 1
			cherrypy.session['perm'] = str(a)
		else:
			cherrypy.session['perm'] = "1"
		cherrypy.session['name'] = "Ales"
		cherrypy.session.save()
		render = web.template.render('/home/ales/integr/website/templates/')
		yield render.tutorial(getStudyTypes(sid).uniqStudyTypes, "")['__body__']
	studyweights.exposed = True

	#####################
	# DISPATCHER function	
	#####################
	def dispatcher(self, **kw):
		cherrypy.session.load()
		sid = cherrypy.session.id
		
		#Make a redirection webpage
		render = web.template.render('/home/ales/integr/website/templates/')
		yield render.testComet()['__body__']
		yield '''<img style="margin-left:auto;margin-right:auto;display:block;" src="07p-6606-loadingcircle.gif" />'''
		yield '''<div id="sliders2">'''
		
		#In case there is session data, the dispatcher will set the session id
		try:
			cherrypy.session.id = kw['sid']
			print (kw['sid'])
			#cherrypy.session.id = kw['sid']
			yield '''Hello, returning user! <br /> Restoring session id... <br /> :) Will redirect you to results <br />'''
			yield 'You will be redirected to results page in '
			time.sleep(1)
			yield '3......   '
			time.sleep(1)
			yield '2......   ' 
			time.sleep(1)
			yield '1......   '
			yield '</div></body></html>'
			cherrypy.session.save()
			yield '''<script>window.location = '/results';</script>'''
		except:
			print ("No SID!")
		
		#In case there is disease choice submitted for static results view, we define this with static=1, if static is not present, an exception is raised and we go to the next step 
		try:
			print (kw['disease'])
			print (kw['static'])
			if kw['static']=="1":
				print ("Haha")
				try:
					tdir = '/home/ales/integr/temp/%s' % sid
					print (tdir)
					if os.path.exists(tdir): 
						try:
							shutil.rmtree(tdir)
							yield "You have some working files in your temporary file storage space.  <br />Cleaning up the temporary dir...  <br /> "
						except: 
							yield "Something wrong with removing the target dir"
					os.mkdir(tdir)
				except:
					print ("There is a problem with defining your temporary storage space!")	
				sdir = '/home/ales/integr/diseases/%s/results.txt' % kw['disease']
				print (sdir)
				if not glob.glob(sdir):
					yield "There is no results for this disease uploaded. Sorry!"
				else: 
					for name in glob.glob(sdir):
						yield "Uploading file %s <br />" % name 
						shutil.copy(name, tdir)
					yield 'You will be redirected to results page in '
					time.sleep(1)
					yield '3......   '
					time.sleep(1)
					yield '2......   ' 
					time.sleep(1)
					yield '1......   '	
					yield '''<script>window.location = '/results';</script>'''
		except:
			print ("No disease name provided!")
		
			
		#In case there is disease choice submitted, the dispatcher will copy the data into current session
		try:
			print (kw['disease'])
			print (kw['static'])
			if kw['static']=="0":
				try:
					tdir = '/home/ales/integr/temp/%s' % sid
					if os.path.exists(tdir): 
						try:
							shutil.rmtree(tdir)
							yield "You have some working files in your temporary dir.  <br />Cleaning up the temporary dir...  <br /> "
						except: 
							yield "Something wrong with removing the target dir"
					os.mkdir(tdir)
				except:
					print ("There is !")	
				sdir = '/home/ales/integr/diseases/%s/Study*.txt' % kw['disease']
				print (sdir)
				if not glob.glob(sdir):
					yield "There is no data for this disease uploaded. Sorry!"
				else: 
					for name in glob.glob(sdir):
						yield "Uploading file %s <br />" % name 
						shutil.copy(name, tdir)
					yield 'You will be redirected to results page in '
					time.sleep(1)
					yield '3......   '
					time.sleep(1)
					yield '2......   ' 
					time.sleep(1)
					yield '1......   '	
					yield '''<script>window.location = '/studyweights';</script>'''
		except:
			print ("No disease name provided!")
			

			
			
		cherrypy.session.save()		
	dispatcher.exposed = True
	dispatcher._cp_config = {'response.stream': True}
		
	def testing(self, **kw):	
		cherrypy.session.load()
		sid = cherrypy.session.id
		print (kw)
		cherrypy.session['studyweights'] = kw
		typeW=[]
		for ustudy in getStudyTypes(sid).uniqStudyTypes:
			print (ustudy)
			typeW.append(int(kw[ustudy]))
			print (typeW)
		cherrypy.session['typeW'] = typeW	
		yield "<br />Parameters, permutation number %s, interval range measuring %s bp as well as study weights are saved for your session id %s. You will be redirected to integration algorithm page!" % (cherrypy.session.get('studyweights')['perm'], cherrypy.session.get('studyweights')['interval'],cherrypy.session.id)
		yield '</div></html>'
		cherrypy.session.save()
	testing.exposed = True

	def test1(self):
		yield "Session id %s" % cherrypy.session.id
		yield cherrypy.session.get('name')
		yield cherrypy.session.get('perm')
	test1.exposed = True
	
	#####################
	# UPLOADER This class will upload the data to the server, using POST
	#####################

	def upload(self, **kw):
		print cherrypy.session.get('name')
		sid = cherrypy.session.id
		
		size = 0
		print (kw)
		while True:
			print (kw)
			fname = kw["Filename"]
			sid = kw["session_id"]
			data = kw["Filedata"].file.read()
			if not data:
				break			
			try:
				tdir = '/home/ales/integr/temp/%s' % sid
				os.mkdir(tdir)
			except:
				print ("Folder already there!")
			tfile = '/home/ales/integr/temp/%s/%s' % (sid, fname)
			if os.path.exists(tfile):
				os.remove(tfile)	
			f = open(tfile, 'w')
			f.write(data)
			f.close()
			return("OK") #Make the scrollbars finish the transfer
		cherrypy.session.id = sid #Correct the session id, that is corrupted by the uploadify flash
		cherrypy.session.save()
	upload.exposed = True
	
	def home(self):
		render = web.template.render('/home/ales/integr/website/templates/')
		yield render.default()['__body__']
		#yield '''<img style="float:right; max-width:35%;" src="Integration_flow-2.png" />'''
		yield '''<div style="float:right;width:240px;height:380px;margin:30px;">''' 

		yield '''<p style="font-size:28px;">Start analyses</p>'''
		yield '''<p> Click on the Start button below to access data on neurodegenerative diseases and to perform custom integrative analyses. <br /></p>'''
		yield '''<input style="margin-left:70px;;margin-top:15px;margin-bottom:15px;background: rgb(41,41,41); " type="submit" class="button" value="Start" onclick="window.location='/index1';"/> '''

		yield '''<p style="font-size:28px;"> Read the guide</p>'''
		yield '''<p> Before beginning the analyses using your custom data we advise reading our manual for using the integratomics tool.  <br /></p>'''
		yield '''<input style="margin-left:70px;margin-top:15px;margin-bottom:15px;background: rgb(41,41,41); " type="submit" class="button" value="Guide" onclick="window.location='/publications';"/> '''

		yield '''</div>'''

		yield '''<div style="margin-top:40px;margin:30px;width:530px;"> <p style="font-size:28px;margin-right:60px;margin-top:20px;">Introduction</p>'''
		yield '''<p style="margin-right:60px;margin-bottom:0px;"> The massive amounts of data, gathered with the development of highly-parallel methods in molecular genetics presented a new challenge to 
				interpret and make new discoveries in these heterogeneous and large biological datasets. We present a new to approach toward tackling this challenge. <br />'''
		yield '''We have created a new framework to integrate genomic data based on genomic location. Results from various sources are projected to genomic intervals and 
				a prioritization approach is used to extract hotspots, where detections accumulate in the genome. We then use permutation approaches to distingish significant 
				accumulation of detection from those resulting from random occurences. <p> '''		
		#yield '''The combination of improving technologies for molecular interrogation of global molecular alterations in human diseases along with increase in computational capacity, have enabled unprecedented insight into disease etiology, pathogenesis and have enabled new possibilities for biomarker development. A large body of data has accumulated over the recent years, with most prominent increase in information originating from genomic, transcriptomic and proteomic profiling levels. The complexity of data, however, made discovery of high-order disease mechanisms involving various biological layers, difficult and therefore require new approaches towards integration of such data into a complete representation of molecular events occurring on cellular level. For this reason, we propose a new mode of integration of results coming from heterogeneous origins, using rank statistics of results from each profiling level. Due to the increased use of next-generation sequencing technology, experimental information is becoming increasingly more associated to sequence information for which reason we have decided to synthesize the heterogeneous results using the information of their genomic position. We therefore propose a novel positional integratomic approach toward studying omic information in human disease. '''

		yield '''<p style="font-size:28px;margin-top:20px;">References</p>'''
		yield '''<p style="font-size:12px;"> Maver A, Peterlin B; Positional integratomic approach in identification of genomic candidate regions for Parkinson's disease. Bioinformatics 2011; 27: 1971-8. <br /></div>'''
	home.exposed=True
	
	def works(self):
		render = web.template.render('/home/ales/integr/website/templates/')
		yield render.default()['__body__']
		
		#yield '''<img style="float:left; max-width:50%;top:300px;position:relative;" src="Integration_flow-2.png" />'''
		yield '''<div style="float:right;width:240px;height:380px;margin:30px;">''' 
		yield '''<p style="font-size:28px;"> Further reading</p>'''
		yield '''<p> For the additional information regarding the algorithm and applications of this algorithm in the field of complex neurologic diseases click below <br /></p>'''
		yield '''<input style="margin-left:45px;margin-top:15px;margin-bottom:15px;background: rgb(41,41,41); " type="submit" class="button" value="Publications" onclick="window.location='/publications';"/> '''

		yield '''<p style="font-size:28px;">Start the analysis</p>'''
		yield '''<p> Click on the Start button below to access data on neurodegenerative diseases and to perform custom integrative analyses. <br /></p>'''
		yield '''<input style="margin-left:70px;;margin-top:15px;margin-bottom:15px;background: rgb(41,41,41); " type="submit" class="button" value="Start" onclick="window.location='/index1';"/> '''
		yield '''</div>'''


		yield '''<div style="margin-top:40px;margin:30px;width:700px;"> <p style="font-size:28px;margin-right:60px;margin-top:20px;">Data integration</p>'''
		yield '''<p style="margin-right:60px;margin-bottom:0px;"> Before data can be integrated, they have to be reduced to a universal common denominator. Due to increasing heterogeneity of genetic information, tying biological information to gene-level annotation is becoming increasingly more difficult. Genomic variation and methylation patterns are two examples of information that is prohibitively difficult to associate with genes in any straightforward manner, as such alterations occur in genes, between genes or spread across several genes. For this reason we opted for an integration based on genomic position of features originating from various data sources. This requires the signals from all datasets to be converted to their genomic positions and projected on the genome assembly backbone. This step then allows for complete omission of difficult annotation conversion steps, required before final integration can be performed, greatly simplifying the synthesis of heterogeneous data. After signals are positioned on the genomic backbone, the complete assembly is divided into bins of equal size. For each study, a score is given to each of the bins, depending on the score of alterations residing in that segment of the genome. After this step, the scores of all bins are prioritized and their rank scores calculated. The integration step is attained when non-parametric rank product for each of the bin is calculated on the basis of rank scores of bins originating from each data source, as we have previously described. The lower final rank product signifies that higher ranks were attained by bins on several separate biological layers. <img style="max-width:100%;" src="Integration_flow-2.png" /> These bins therefore represent genomic regions where accumulation of signals is detected on various biological levels and thus represent regions of interest for further investigation. Ultimately, a permutational test may be employed to determine the significance of signal accumulation in each bin. The detailed overview of the process may be observed in Figure 1.'''
		yield '''We have created a new framework to integrate genomic data based on genomic location. Results from various sources are projected to genomic intervals and 
				a prioritization approach is used to extract hotspots, where detections accumulate in the genome. We then use permutation approaches to distingish significant 
				accumulation of detection from those resulting from random occurences. <p> '''		
		#yield '''The combination of improving technologies for molecular interrogation of global molecular alterations in human diseases along with increase in computational capacity, have enabled unprecedented insight into disease etiology, pathogenesis and have enabled new possibilities for biomarker development. A large body of data has accumulated over the recent years, with most prominent increase in information originating from genomic, transcriptomic and proteomic profiling levels. The complexity of data, however, made discovery of high-order disease mechanisms involving various biological layers, difficult and therefore require new approaches towards integration of such data into a complete representation of molecular events occurring on cellular level. For this reason, we propose a new mode of integration of results coming from heterogeneous origins, using rank statistics of results from each profiling level. Due to the increased use of next-generation sequencing technology, experimental information is becoming increasingly more associated to sequence information for which reason we have decided to synthesize the heterogeneous results using the information of their genomic position. We therefore propose a novel positional integratomic approach toward studying omic information in human disease. '''

		yield '''<p style="font-size:28px;margin-top:20px;">References</p>'''
		yield '''<p style="font-size:12px;"> Maver A, Peterlin B; Positional integratomic approach in identification of genomic candidate regions for Parkinson's disease. Bioinformatics 2011; 27: 1971-8. <br /></div>'''
	works.exposed=True
	
cherrypy.root = HelloWorld()

current_dir = os.path.dirname(os.path.abspath(__file__))

cherrypy.config.update({
	'log.screen':True,
	'tools.sessions.on': True,
	'tools.sessions.storage_type': 'ram',
	'tools.sessions.timeout': float(360),
	#'tools.sessions.locking' = 'explicit',
	'checker.on':False,
	'server.socket_host': '0.0.0.0',
	'server.socket_port': 21,
	})
	
conf = {
	'/': {
		'tools.staticdir.on' : True,
		'tools.staticdir.root': '/home/ales/integr/website/static/uploadify',
		'tools.staticdir.dir' : ''},
	'/index1':{
		'tools.staticdir.on' : True,
		'tools.staticdir.dir' : ''}}


#cherrypy.tree.mount(HelloWorld())
cherrypy.tree.mount(HelloWorld(), config=conf)
cherrypy.engine.start()
cherrypy.engine.block()

