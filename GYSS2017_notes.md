#GYSS 2017 summit notes
---

###Talk 1 - Online data and personal control - Butler Lampson
what makes online data special:

- widespread
- accessible
- data on the physical world

One challenge is that technology is tied to governance and rules

- different governments have different rules

data and data tags (metadata) that defines the rules of data use policy

__rule 1 - data must have policy copied as well__

policy service online

agent / handeler / policy service online

this whole thing only works for agesnts that are subject to regulation, in which there are many playners that are not.

```If you can't do anything about rogue users, what the hell is the point```

It also does not protect perople from coercion, which is a serious problem.

How to work across all devices, services, and anything else that is available...

```This is expecting that people care```

Who controls the policy services, they have the power.

At what level can you remove data policy tags when data is aggregated?

---

###Talk 2 - can we cure all diseases soon? - Aaron Ciechanover

diseases easier to prevent than to treat.

Various drug revolutions

__first rev__ - the era of incidential discoveries (1930-1960s)

- asprin
- penicillin

__second rev__ - high throughput bute force screening of large libraries of chemical compounds

- robot based chemical screens
- statins

__21st century__ - personalized, predictive, preventive, participatory medicine (the 4P's of medicine)

- target medicine and treatment protocols that are taylored to you specifically, as we are not all the same
- stem cell-based therapies
- Herceptin (targeted in breast cancer)

---

###Talk 3 - The power of abstraction - Barbara Liskov

What lead us to structuring software the way we do today?

Software today is very complex

To address complexity, we have developed a ariety of algorithms, data structures, and protocols. In the beginning, none of these existed.

The design and structure of the software requies to fit your needs, this is still hard and a requirement even if you use many of these concepts that exist (__programming methodology__ leading to programming langauges)

In 1970's - the software crisis

_Dijkstra Cacm 1968 - Go To Statement Considered Harmful_

"reasoning dynamic (program output) from static representation (code) is very hard when using goto statements"

_Parnas IFIP 1971 - Information Distribution Aspects of Design Methodology_

"The connextions between modules are the assumptions which the other modules make about each other"

Programs are collections of modules, each of which has an __interface__ described by a __specification__ plus __implementations__

Modularity provides 

- local reasoning
- modifiability (you can throw one out and replace them)
- independent development (with standardized specificaitons, things can be built seperatly and be expected to come together in the end)

Big idea: Connecting partitions to data types

Partition state to represent an object which is connected to opperations that are important to that object data type

What defines data types?

_Dahl & Hoare Structured Programming 1972 - Hierarchical Program Structures_

_Morris Cacm 1973 - Protection in Programming Languages_

"no program on the otuside can modify the state internally"

"no program on the outside should even look at the internal state"

How can you actually make this happen?

Abstract data dypes - a sketch of a programming language (CLU)

- Clusters
- Polymorphisms, e.g. set[T]
- Iteration
- Exceptions

Followed by Argus and distributed computing

Modularity based on abstraction is the way things are done in programming today (it is just "the way it is") - hence the value of her work!

---

###Talk 4 - A computation Theory of the Corext - Leslie Valiant

__Turing's Thesis (1936)__ - What is mechanistically computable (in any reasponable intuitive sense) is exactly what is computable by a Turing Machine.

If the brain, evolution, biology are not computations in this sense _then Turing's thesis is wrong._

More recent studies of mechanisms also measure resourees or complexity, such as the number of steps.

To apply this: Any suggested machanism for a phenomenon is wrong if it needs more steps than the phenomenon takes.

A first model of the brain as a neural net (1943) - McCulloch & Pitts

Neuroscience doesn't like top down theory-driven work. They like observations!

How does the brain work? - we don't know.

Symbolic processing

The cortex can do this for hundres of thousands of such taslks in succession - but we dont' know how it can even in principle

So what do we do wheen there are no viable theories to test?

1. Use neural model that underestimates the brain
2. Specify some challenging set of multiple task types
3. Show that these task types can be executed on the model
4. Show that sequences of thousands of interleaved tasks of these types can be supported without degradation of earlier ones.

__Is there anything simple about the brain?__

Yep - the brain is slow and communication challenged (each neuron is connected to only a small faction of others). No global addressing in neurons. The effect of neuron on others is quite weak.

How many neurons are required to influence another to fire

k=1 single leader
k=lots = super democratic decision

larger k makes things more difficult

How does the brain allocate storage to a new concept?

__Hierarchical Memorization__: For any stored items A, B, allocate neurons to new item C and change  synaptic weights so that in future A and B active will cause C to be active also. (Chunking)

This has been shown in part in experiments...

How do we relate and make relationships among represented concepts?

__Association__: For any stored items A, B, change synaptic weights so that in future when A is active then B will be caused to be also.

Can create models:

Nodes (n), edges (d) , weights (1/k) x threshold

```started to get a bit deep here...```

Can begin to perform experiments to test these models (Connectomics, Optical stimulation and recording, MEG)

---

###Talk 5 - Reconstruction of the Past - Sydney Brenner

One of the biggest struggles of science today is understanding what has happened in the past.

From theories of Natural selection to our understaindg of DNA and heritable information, we have gained greater understanding of the evolution process. Soon, we can create a new synthesis that ties all life together.

random selection vs survival of the fittest

what brought chimpanzees into humans? A big question for todays scientists

man survives due to having an inventive brain?

---

#Day II

###Talk 6 - What is constructive mathematics and why it is important - Vladimir Voevodsky

Invesnt or learn mathematican constructions andthen use them to prove theorems that answer questions that has been posed by other, mostly seniorl mathematicians.

Math = 

- definitions
- constrictions
- assertions about these constructions
- and proofs of these assertions

As the field is based on building on history, errors historically can lead to 'crashes' in the field that break down a wide variety of gains made over these flawed solutions.

__consrtuctive vs non-constructive assertions__

non-constructive = 'classical'

```constructive = proofs?```

What is gained by constructives and proofs?

math is deeply tied in an abusive relationship with religion

_"The belief in methematics as the basis of the structire of our world. From a simple one of Pythagoreans to string theory."_

The importance of mathematical reality as the ultimate reality on which the everyday reality is built.

__Different methods of proof and construction that are permitted by classical and constructive mathematics lead to different mathematical realities.__ Which one of these is the ultimate one behild out world?

Math provides the tools and language to define our reality 

---

###Talk 7 - Two Pillars of Physics: Quantum Mechanics (90 years ago) & General Relativity (100 years ago) - David Gross

The birth of quantum mechanics

planke - einstein - bohr

followed by revolution in final formulation (heisenberg)

The state of quantum mechanics is that

- It works
- It makes sense
- it is hard to modify

So that's good.

replace phase space with 'hilbert space' (vectors and rays)

In quantum mechanics one can describe a system in may incompatible different ways. __There is no unique exhaustive description__

__Will it eventually fail somewhere?__

The 'collapse of the wave function' upon measurement occurs within the human mind.

- large complex systems?
- conscious systems?
- small distances?
- in describing the whole universe?

so far, it works from atoms, molecules, to large bodies of things. Some examples where it works:

- superconductivity (1911)
- spontanious breaking of E & M (1957)
- ElectroWeak Symmety Breaking
- Exotic metals, heavy elections, spin liquids, etc

Also works at the __subatomic level__

Quantum field theory - you get extreme precision in measurements, but lacking in rigor. But still leads to accurate results

Ranges from 10^27cm to 10^-18cm in which experiments hold true. Expect it to work down to 10^-33cm

3 families of quarks and leptions

e & M / Strong / and weak force

across 60 orders of magnitude it has been shown. ```neat```

It basically dominates modern technology to explain our world:

- transistor
- resistor
- laser
- integrated circuit
- MRI / PET / other medical imaging

__What are the future prospects?__

- the discovery and or creation of new forms of matter (these are exciting)
- influence at the atomic level (graphene, nanotubes)
- quantum computers (moving from ones and zeros and 'up' and 'down' into non-binary states

The environment leads us to classical reality.

__General relativity__

mass tells space-time how to curve, and space-time tells mass how to move

```Gen. relativity = Newtonian Gravity + Gravitational Waves```

- dynamic space-time
- physical cosmology (we have a good model from start to now in universe)
- unified theory (still undetermined)

unifying quantum and relative have some issues.

Forces unify where gravidy reach a strong force

string theory (are all aprtiles different vibrations of a superstring?)

unifying quantum mechanics and spacetime is big goal

current view of spacetime

- smooth manifold
- fixed topology
- fixed number of dimensions

but is falling apart!

Is all of spacetime emergent?

Is gravity emergent?

_what is spacetime made of?_

So the last 100 years still drives fundimental physics and the best is yet to come!

---

###Talk 8 - Towards Next generation Antibiotics - Ada Yonath

Let's talk Ribosomes, found in all cells. ribosomes act fast and are accurate

- 40 peptide bonds in a second

small and large subunit that holds mRNA and binds in tRNA molecules

peptide bonding in large subunit and built in tunnel for protection. Once outside the unit and tunnel, it can begin folding accordingly.

as ribosmes are so fundimental, many antibiotics target them

over 40% of clinically useful antibiotics target protein biosynthesis, mostly by paralyzing the ribosome.

Natural antibiotics are weapons that bacteria from one type use for interfering with cell lofe of different species of microorganism

How do small antibiotics paralzse the giant ribosome (500-100Da vs 20k+)

They bind to functional sites in both large and small subunits.

- tunneling
- decoding
- hinging
- bonding

__A big problem__ is that functional sites are targets bit are highly conserved across species. How can use make sure they only are useful on pathogenic bacteria and not our own ribosomes...

Must define by the subtle differences that exist as ribosomal gene mutations

__Antibiotic resistance__ can be created by mutations that influence those functional sites. It is hard to prevent this from happening.

A need to move to understanding pathogen structures to strain-specific structure

Outcomes:

1. The eficency of antibiotics can often be imprived
2. We can imprive antibiotic selectivity

ID non-functional sites that influence functional sites (look for indels, dna changes in ribosomes) can design antibiotics to impact these targets which have not bee nselected for resistnace so far.

Most ribosomal antibiotics are extensions of small organic molecules that cannot be digested by eukaryotes and cannot be degradable in the environment.

Can design new antibiotics that are environmentally degradable.


"Pathogen specific antibiotics'

ID structures quickly via 3d cryo electron microscopy compared to old x-ray crystallography

goals are to __minimize__ wide spread resistance.

---

###Talk 9 - What can we do with a quantum liquid? - Anthony Leggett

super cold physics

Low temps provide a high degree of order

Quant Mechanics provide particles that show 'wave' behavior

Need liquis that allow particles to change places

Quantim liquids show effects of quantum mechanics plus indistinguishibility

thesse lead to superfluids

###Talk 10 - Energy Beyond Old - Muchael Gratzel

10kw per person in US for lifestyle

---
---
---
###Talk 11 - The Internet of Things - Vint Cerf

"Cyber-physical systems"

Any piece of software that takes real world input, does something, and then has effect on the world.

How far can this concept be expanded?

Web-langauges?, any programable device??

__Biggest concerns:__

- Bugs (who will fix? __for how long?__)
- Updates (source? integrity? confirmations?)
- (Strong) authentication of parties (who can control? configure? receive data?)
- Scaling up configuration (how to do this automatically in secure ways)

Machine learning and AI

- __"The artificial idiot"__ that learns the wrong things
- Ephemeral or Episodic Access (Fire, police, medical responses)

Access control and security

- Residents vs Guest rights (residents vs guests vs robbers?, what happens when they leave?)
- Parents vs Kids
- Show floor operations vs Administrative staff and visitors
- Access to data gathered by device in operation

The role of standards

- Major issue of interoperability, safety, privacy, and security
- Proprietary and private sector consortia conventions
- - IEEE, THREAD, WEAVE, SGIP, ZIGBEE, IOTX, Schema.org, OCF, IPv6, NEMA, etc..
- National and International Standards (ANSU, ITU-T, ISO...)
- Backward compatibility
- Likely pressure from consumers, business, government, miliraty, health care, ...

```will there be a backlash by noting the benifits of non-IOT things??```

Long term questions

- Digital picture frames (storage corruption and power supply failure)
- Custom systems (entertainment, lighting) - bigs and other protocol dependencies
- Programmable light switches (battery format dependency, programing UI is bad)

the dependency chain of items required to work...

Major considerations

- reliability and ease of use
- safety
- security
- privacy
- autonomy (what can work when the internet is down?)
- interoperability (of ensembles)

cyber underwriters lab

---

###Talk 12 - Switches & Latches: The control of cell Division - Sir Timothy Hunt

how does one cell become two cells

```The biology of the cell cycle - J.M. Mitchison```

replication and segregation

M - G1 - S - G2 - phases

Three switches, points of no return

 - Enter mitosis (all DNA is replicated)
 - Get out of mitosis (chromosomes are all ligned up)
 - starting DNA replication

 Research on these 'switches' by looking at Urchin Eggs
 
 protein synthesis increases after fertilization
 
 Maturation Promoting Factor (MPF) enzyme that catalyzes cell division. Found in all cells and highly conserved!
 
 something turns MPF on and off
 
 __cyclins__ go up and down in patterns related to cell division, could they do the trick?

cdc2 and cdc28 

 - controls G1 to S and G2 to M
 - encodes protein kinase that is inactive
 - seems to be tied to cyclins as well (required for activated version!)
 - CDC2 + CYCLIN = MPF

Protein Kinase activates from phosphorilation leading to 100x more activity

cell makes cyclin, cdc2 turns on, enters mitosis

cyclin is degraded, cdc2 turns off, exits mitosis

hundreds of targetes (or more) that MDF phosphorylates, some targets get multiple P

```The futile cycle problem```
kinase adds, phosphtase removes, can't use both at the same time to get anywhere

phosphotase turns off in mitosis by being inhibited by a pseudo-substrate

two positive feedback loops

CDK kinase promotes Gwl (great wall) which phosplorylates the inhibitor that stops the phosphotase working. This moves out of mitosis

---

###Talk 13 - The traveling-salesman problem - Richard Karp

minimize total travel cost between cities given known cities and distances between

searching combintorial possibilities to find the best one. A genreal problem

can be viewed as symmetric costs between nodes ```Cij = Cji``` or asymmetric 

Can be solved for completion mathematically, have been done for 80k+ cities

__The problem is NP-complete__

This problem is solvable in polynomial time if and only if ```P=NP``` where ```P``` is the class of problems solvable in polynomial time, and ```NP``` is the class of roblems for which solutions can be verified in polynomial time.

if we don't have a polynomial time algorith for this, can we use alternatives to obtain an exact solution?

related problems

- Minimum-Weight Spanning Tree (find minimim weight of edges to connect all)
- Minimum-Weight Matching (pair up cities to find the minimum weights)
- Eulerian Graph (a connected graph with all degrees even) - can convert any Eulerian weighted graph to a tour whoes total weight is less than the sum of the weights in the original graph.

Christofides Approximation Algorithm can approximate the traveling salesman problem using these above techniques

The total cost of the resulting tour is at most 3/2 times the cost of an optimal tour (haven't found anything less than this to date)

As far as empirical results (practical value), there are algorithms that can prodice a tour cost at most 1.02 times the cost of an optimal tour by looking for local improvements.

Linear programming is required to truly find the optimal tour.

`so lots of math to solve math things`

---

###Talk 14 - Random Walk to Graphene - Sir Andre Geim

How and why

random experiments to get there inthe end

magnetic water descaler - what do magnets do to water?

levetating water in magnetic field - indicates the ever present diamagnetism which is generally small, but not as negligible as commonly believed

lots of things in the magnetic field (including frogs), lots of things levetate

Why are geckos able to climb up walls? - sticky feet with lots of interactions and hairs at the micron scale. Can hold things up quite greatly

Another ranodm idea: Make films of graphite as thin as possible and study their properties

Materials indicate that nature hates low dimentions - growth of one atom thick materials is forbidden in nature...

...but does not mean it cannot be made artificially.

take 3d materia and pull out individual atomic plane

__Graphene__ is an extremely simple structure

many unique properties - many superlatives to describe it

very special as you can mimic ultra-relativistic particles which is generally only seen near the speed of light.

Klein Tunnelling

Impacts the fundamental limits of the periodic table by pushing beyond the 116th element (generally limited by the supercritical regime)

- mimics relativistic physics to create 'artificial atoms' which can overcome the overcriticle state

So lots of unusual effects are known and more are being found

__What are its applications?__

lots

Other 2d crystals have been made, so many different are being tried and examined for unique properties.

make stacks of different 2d crystals to make materials on demand

lots of research here

"New tool in our human toolbox"

some flakes from your graphite pencil are going to be graphene - we've always had it around is. `highlights how little we still know about the world around us`

---
###Talk 15 - Dairy Cattle Products as Risk Factors for Human Cancers - Harald zur Hausen

~21% of global cancer incidence is linked to infections

Lots of cancers have been sequenced, however no consistent reports appeared pointing to novel infections agents in human tumors.

The concequence is that we have either have discovered the majority of infections, __or__ the involvement of infectious agents does not follow the pattern of well established oncogenic pathogens

High risk regons for colon cancer seem to be linked to dairy fbeef concumption. Pork, poultry, fish diets seem to play no significant role.

`stop eating food`

species-specific factor in dairy beef (specific bacteria have been reported to be involed in colon polups, cancer, and other inflammatory reactions

cow macteria and viruses can be passed to people. This is bad.

over lifetime, consumption of beef will induce mutations in host cell DNA, this mutational load can lead to  colon cancer. (Transient or latent infection by carcinogenic virus

may be tied to other cancers as well.

`this is depressing`

`prob should go vego`

breast cancer and bovine milk factor

prolonged breastfeeding has protective effect (>6mo) on chronic diseases

---

###Talk 16 - The Challenge of Reconciling the Physics of Black Holes with Quantum Mechanics - Gerard't Hooft

---

###Talk 17 - Perspectives on Mathematics and Biology - Stephen Smale

1950's = botany and zoology was divide in the field

Now, math is a big part

- computational biology
- machine learning

talked about chromatin a bit

"fixed network with dynamic chromatin layer"
---
###Talk 18 - Beyond charge currents: spon and ion currents for future computing technologies - Stuart Parkin

Moving beyond the electron for computing!

in last 100 years, we have used a variety of technologies for computing:

- mechanical
- electro-mechanical
- vacuumtube
- transistor
- integrated circuit

we are coming to the end of silicon technology and moving to:

- quantum computing (restrictive usecases)
- cognitive computing

moving from charge to orbitals, spin, ions

__spintronics__ - from materials and phenomena to applications

hd's are slow
MRAM - development based on magnetetic states - a 2D technology

__Magnetic Racetrack Memory__ - currents to move magnetic sections across sensors

- 3d storage class
- solid state
- capacities of current HDs x 1mill faster
- can be expanded into spin polarized current to manipulate domain walls

`I'm retarded`

magnetic domain-wall nanowire shift register has been shown to work

