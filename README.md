# digital twin simulating uterine myomas/fibroids
uterine fibroids affect 70% of women before they reach menopause, with around 40% of women being symptomatic, meaning they experience significant pain, heavy bleeding, and fertility issues. most treatments to this issue include hormone therapy and invasive surgeries, and these hormone therapies can cause numerous detrimental side effects. despite these challenges, the research done in the women's healthcare field is lacking and especially doesn't explore more natural and less invasive treatments to these reproductive issues. this project is an attempt to simulate how potential treatments that are less researched impact the fibroids.

## data sources 
the MRI dataset comes from the UMD dataset, linked [here](https://www.nature.com/articles/s41597-024-03170-x). i was not, however, able to find good datasets showing sex hormones across a menstrual cycle and certainly wasn't able to find a categorical dataset of any sort with significant or helpful information for this project, in part due to a lack of study or due to data sharing restrictions :(. to compensate for this, i've built the hormonal dataset and categorical treatment dataset synthetically and based upon published physiology. i injected some sparsity/missigness within the data so it's realistic and emulates real, questionably reliable data, so some data points that would be self-reported or retrieved from blood tests are missing.

## approach 
this project has two segments, both built upon largely synthetic/generated datasets:

#### initial multimodal framework for building a predictive simulation to model tumor growth
i have three datasets, and to represent the full biological environment, i'm training a model on each dataset separately and then merging the results:
- hormonal time series dataset: tracks womens' hormone levels across each phase of her cycle
  - this is important because estrogen and progesterone levels across a woman's cycle directly impact fibroid growth
  - this dataset was synthetically generated. more about this in the following sections
- mri imaging dataset: includes 3d mri / voxels of an mri scan, labeled with fibroid tissue vs healthy tissue
  - this data describes real fibroids and is the only real source of truth
- categorical treatments dataset: represents women with varying severity of fibroids and the corresponding treatment plan that they are currently on
  - bridges the two by linking patient profiles, fibroid severity, and treatment plans together
  - this dataset was synthetically generated. more about this in the following sections
  
the three datasets share overlapping feature spaces (fibroid severity, treatment type, patient demographics) that enable cross-modal reasoning. the categorical dataset bridges MRI-derived volume measurements with treatment plans, providing the context needed for the PINN to model growth trajectories under different interventions.

generally, this project is intended as a proof of concept for what a multi-modal fibroid environment model would look like with sufficient real data.

#### predicting tumor response to treatment using physics-informed neural networks
i synthesized and generated a dataset that shows tumor size/volume over time, and its growth or shrinkage is dependent on how aggressive the treatment applied is. the physics-informed neural network (pinn) on top of this dataset extrapolates conclusions and makes educated guesses on tumor growth or shrinkage based on treatment.

the verhulst logistic growth model equation is:

$$
\frac{dV}{dt} = rV(1 - \frac{V}{K})
$$

where r = growth rate, and K = carrying capacity

modifying it so it can also model tumor decay/shrinkage:

$$
\frac{dV}{dt} = rV(1 - \frac{V}{K}) - kV
$$

where k = decay rate.

i used the closed form of the decay function to generate the synthetic dataset of tumor growth/shrinkage based on different treatments, which have different corresponding decay rates. the residual of this equation is used when finding the loss function:

$$
loss = MSE(V) + \lambda MSE(residual),   \text{  where residual} = \frac{dV}{dt} - rV(1 - \frac{V}{K}) + kV
$$

## preprocessing 
the MRI volumes are downsampled to 128×128×5, which is small enough to keep compute reasonable, but large enough to preserve clinically meaningful fibroid geometry without losing smaller lesions. t2 volumes use bilinear interpolation to maintain smooth intensity gradients, while segmentation masks use nearest-neighbor interpolation to preserve exact label boundaries, since blending between "fibroid" and "healthy" tissue labels would produce meaningless intermediate values. t2 volumes are z-score normalized per volume to account for scanner intensity variation across different scans. fibroid presence, count, and volume ratio relative to total uterine tissue are extracted from the segmentation masks using connected component analysis!

the hormonal and categorical datasets are normalized using standard scaling and one-hot encoded for categorical variables.

`preprocess.py` implements this and also generates the datasets.

## models
`mlp_categorical.py` is a simple binary classification model that determines whether someone likely has fibroids or not, based on a categorical dataset containing data for pain level, cycle length, ferritin level, etc. this dataset is intentionally sparse so it can be used for an ml systems project that decides when the model should abstain from making a decision at inference time, due to data unreliability. project is linked [here](https://github.com/vanigupta123/data-decision-maker).

this model has an accuracy of 96.61%.

`pinn.py` implements the pinn described above in approach > "predicting tumor response to treatment using physics-informed neural networks".

this model has the following results: mse=0.0266, mae=0.1129, r2=0.1911.
## limitations & honest framing 
while this project explores how potential treatments may impact fibroids, it doesn't take any of the woman's other body systems into account.

as mentioned above, two of the three datasets are synthetically built from published physiology rather than real patient data. this means the models trained on hormonal and categorical data will reflect the assumptions baked into the synthesis process, so they can't discover relationships that weren't already encoded during data generation. the MRI component is the only one grounded in real patient data and carries the most scientific weight. treatment impact predictions are not currently meaningful and shouldn't be interpreted as such. the long-term goal is to incorporate real hormonal and treatment outcome data as it becomes available and to extrapolate behaviors using physics-informed constraints based on known endocrine dynamics.

## setup/usage 
```
git clone https://github.com/vanigupta123/uf_digital_twin
cd uf-digital-twin
pip install -r requirements.txt
python preprocessing/preprocess.py
```
then run whichever model you want within `\models`! ex:
```
python models/pinn.py
```
note: some of this will break locally unless you download the UMD dataset. info for that is in the `data sources` section above.
