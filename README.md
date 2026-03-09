# digital twin simulating uterine myomas/fibroids
uterine fibroids affect 70% of women before they reach menopause, and they cause significant pain, heavy bleeding, and pose a threat to fertility. most treatments to this issue include hormone therapy and invasive surgeries, and these hormone therapies can cause numerous detrimental side effects. despite these challenges, the research done in the women's healthcare field is lacking and especially doesn't explore more natural and less invasive treatments to these reproductive issues. this project is an attempt to simulate how potential treatments that are less researched impact the fibroids.

# approach 
i have three datasets, and to represent the full biological environment, i'm training a model on each dataset separately and then merging the results:
- hormonal time series dataset: tracks womens' hormone levels across each phase of her cycle
  - this is important because estrogen and progesterone levels across a woman's cycle directly impact fibroid growth
  - this dataset was synthetically generated. more about this in the following sections
- mri imaging dataset: includes 3d mri / voxels of an mri scan, labeled with fibroid tissue vs healthy tissue
  - this data describes real fibroids and is the only real source of truth
- categorical treatments dataset: represents women with varying severity of fibroids and the corresponding treatment plan that they are currently on
  - bridges the two by linking patient profiles, fibroid severity, and treatment plans together
  - this dataset was synthetically generated. more about this in the following sections
  
using all three datasets creates a shared context, which is necessary when merging the three models.

while i complete working on the models, this project still serves as a complete preprocessing pipeline! generally, this project is intended as a proof of concept for what a multi-modal fibroid environment model would look like with sufficient real data.

# data sources 
the MRI dataset comes from the UMD dataset, linked [here](https://www.nature.com/articles/s41597-024-03170-x). i was not, however, able to find good datasets showing sex hormones across a menstrual cycle and certainly wasn't able to find a categorical dataset of any sort with significant or helpful information for this project, in part due to a lack of study or due to data sharing restrictions :(. to compensate for this, i've built the hormonal dataset and categorical treatment dataset synthetically and based upon published physiology. 

# preprocessing 
the MRI volumes are downsampled to 128×128×5, which is small enough to keep compute reasonable, but large enough to preserve clinically meaningful fibroid geometry without losing smaller lesions. t2 volumes use bilinear interpolation to maintain smooth intensity gradients, while segmentation masks use nearest-neighbor interpolation to preserve exact label boundaries, since blending between "fibroid" and "healthy" tissue labels would produce meaningless intermediate values. t2 volumes are z-score normalized per volume to account for scanner intensity variation across different scans. fibroid presence, count, and volume ratio relative to total uterine tissue are extracted from the segmentation masks using connected component analysis!

the hormonal and categorical datasets are normalized using standard scaling and one-hot encoded for categorical variables.

# limitations & honest framing 
while this project explores how potential treatments may impact fibroids, it doesn't take any of the woman's other body systems into account.

as mentioned above, two of the three datasets are synthetically built from published physiology rather than real patient data. this means the models trained on hormonal and categorical data will reflect the assumptions baked into the synthesis process, so they can't discover relationships that weren't already encoded during data generation. the MRI component is the only one grounded in real patient data and carries the most scientific weight. treatment impact predictions are not currently meaningful and shouldn't be interpreted as such. the long-term goal is to incorporate real hormonal and treatment outcome data as it becomes available and to extrapolate behaviors using physics-informed constraints based on known endocrine dynamics.

# setup/usage 
```
git clone https://github.com/vanigupta123/uf_digital_twin
cd uf-digital-twin
pip install -r requirements.txt
python preprocessing/preprocess.py
