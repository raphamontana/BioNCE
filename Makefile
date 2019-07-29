data:
	python src/data/make_dataset.py
	python src/features/build_features

train:
	python src/models/train_model.py
