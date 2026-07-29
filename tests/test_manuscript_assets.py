from mapt_zero_shot.manuscript_assets import pearson


def test_pearson_perfect_positive_and_negative():
    assert abs(pearson([1, 2, 3], [2, 4, 6]) - 1.0) < 1e-12
    assert abs(pearson([1, 2, 3], [6, 4, 2]) + 1.0) < 1e-12
