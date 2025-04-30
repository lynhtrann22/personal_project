### Background / Context:
The WiDS Datathon Global Challenge was developed with Ann S. Bowers Women’s Brain Health Initiative (WBHI) in collaboration with Cornell University and UC Santa Barbara. Datasets and support are provided by the Healthy Brain Network (HBN), the signature scientific initiative of the Child Mind Institute, and the Reproducible Brain Charts project (RBC).

Neuropsychiatric disorders that occur in development, like anxiety, depression, autism, and attention deficit hyperactivity disorder, or ADHD, often differ in how and to what extent they affect males and females. ADHD occurs in about 11% of adolescents, with around 14% of boys and 8% of girls having a diagnosis. There is some evidence that girls with ADHD can often go undiagnosed, as they tend to have more inattentive symptoms which are harder to detect. Girls with ADHD who are undiagnosed will continue suffering with symptoms that burden their mental health and capacity to function.

### Challenge Overview:
In this year’s WiDS Datathon, participants will be tasked with building a model to predict both an individual’s sex and their ADHD diagnosis using functional brain imaging data of children and adolescents and their socio-demographic, emotions, and parenting information.

#### Challenge Question and Task:

“What brain activity patterns are associated with ADHD; are they different between males and females, and, if so, how?”

The task is to create a multi-outcome model to predict two separate target variables: 1) ADHD (1=yes or 0=no) and 2) female (1=yes or 0=no).

### Evaluation
The F1 score is the harmonic mean of the precision and recall. It thus symmetrically represents both precision and recall in one metric. The highest possible value of an F-score is 1.0, indicating perfect precision and recall, and the lowest possible value is 0, if precision and recall are zero.

Since the theme of this challenge is to uncover gender inequities and because ADHD diagnosis is harder for females to predict, for the purposes of this competition we are assigning 2x weight to Female ADHD cases (ADHD_Outcome=1, Sex_F=1). In our implementation of the F1 Score, weighted F1 Score is calculated on each column, and those two individual scores are averaged to get the final Kaggle leaderboard score.
