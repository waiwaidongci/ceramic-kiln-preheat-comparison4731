"""瓷砖预热独立批次的 Welch 检验与质量放行护栏。"""

import math
import statistics
from datetime import date
from fractions import Fraction


STANDARD = "现行"
STAGED = "分段"
RECIPES = (STANDARD, STAGED)
REVIEW_PENDING = "待核对"
REVIEW_CONFIRMED = "已核对"
REVIEW_DEPENDENT = "发现配对或相关批次"


def _text(value, label):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label}不能为空")
    return value.strip()


def _number(value, label, minimum, maximum=None, *, whole=False):
    message = f"{label}须为整数" if whole else f"{label}须为有限数字"
    if isinstance(value, bool) or not isinstance(value, (int, float, str)):
        raise ValueError(message)
    try:
        number = int(value) if whole else float(value)
    except (ValueError, OverflowError) as error:
        raise ValueError(message) from error
    if (whole and isinstance(value, float) and number != value) or (
        not whole and not math.isfinite(number)
    ):
        raise ValueError(message)
    if number < minimum or (maximum is not None and number > maximum):
        if maximum is None:
            raise ValueError(f"{label}不能低于 {minimum}")
        raise ValueError(f"{label}须在 {minimum} 到 {maximum} 之间")
    return number


def _parse_batches(rows):
    groups = {recipe: [] for recipe in RECIPES}
    batch_ids = set()
    scope = None
    for index, row in enumerate(rows, start=1):
        prefix = f"第 {index} 条批次记录"
        batch_id = _text(row.get("batch_id"), f"{prefix}的生产批号")
        if batch_id.casefold() in batch_ids:
            raise ValueError(f"生产批号重复：{batch_id}。每个独立批只能录入一次。")
        batch_ids.add(batch_id.casefold())
        prefix = f"生产批 {batch_id}"
        raw_date = _text(row.get("run_date"), f"{prefix}的生产日期")
        try:
            run_date = date.fromisoformat(raw_date).isoformat()
        except ValueError:
            raise ValueError(f"{prefix}的生产日期须为有效的 YYYY-MM-DD") from None
        if run_date != raw_date:
            raise ValueError(f"{prefix}的生产日期须为有效的 YYYY-MM-DD")
        kiln = _text(row.get("kiln_id"), f"{prefix}的窑线")
        product = _text(row.get("product_code"), f"{prefix}的砖种规格代码")
        if scope is None:
            scope = (kiln, product)
        if scope != (kiln, product):
            raise ValueError(f"{prefix}超出对照范围。请只导入同一窑线、同一砖种规格的批次。")
        recipe = _text(row.get("recipe"), f"{prefix}的预热工艺")
        if recipe not in groups:
            raise ValueError(f"{prefix}的预热工艺须填写「现行」或「分段」")
        area = _number(row.get("area_m2"), f"{prefix}的投入面积", 0)
        energy = _number(row.get("preheat_kwh"), f"{prefix}的预热段电量", 0)
        if not area or not energy:
            raise ValueError(f"{prefix}的投入面积和预热段电量均须大于 0")
        intensity = energy / area * 1000
        if not math.isfinite(intensity) or intensity == 0:
            raise ValueError(f"{prefix}的电耗强度超出可计算范围，请核对电量和面积的单位")
        inspected = _number(row.get("inspected"), f"{prefix}的抽检片数", 1, whole=True)
        defective = _number(row.get("defective"), f"{prefix}的缺陷片数", 0, inspected, whole=True)
        groups[recipe].append((run_date, batch_id, intensity, inspected, defective))

    if not all(groups.values()):
        raise ValueError("台账须同时包含现行和分段预热工艺的批次")
    return scope, groups


def validate_rows(rows):
    """Check imported data before it is saved to the batch ledger."""
    _parse_batches(rows)


def _get_settings(params):
    alpha = _number(params.get("alpha", 0.05), "单侧显著性水平 α", 0.001, 0.2)
    saving = _number(params.get("min_saving", 20), "最低节能量", 0)
    min_batches = _number(params.get("min_batches", 8), "每组最少独立批数", 3, whole=True)
    min_inspected = _number(params.get("min_inspected", 1200), "每组最少抽检片数", 1, whole=True)
    max_rate = _number(params.get("max_defect_rate", 1.5), "分段组缺陷率上限", 0, 100)
    max_increase = _number(params.get("max_defect_increase", 0.3), "缺陷率涨幅上限", 0, 100)
    review = params.get("design_review", REVIEW_PENDING)
    if review not in (REVIEW_PENDING, REVIEW_CONFIRMED, REVIEW_DEPENDENT):
        raise ValueError("请选择有效的独立批次试验核对状态")
    return alpha, saving, min_batches, min_inspected, Fraction(str(max_rate)), Fraction(str(max_increase)), review


def _regularized_beta_lower(a, b, x):
    # Modified Lentz evaluation of the regularized incomplete beta expression (DLMF 8.17.22).
    floor = 1e-300
    previous_denominator = 0.0
    previous_numerator = 1.0
    denominator_fraction = 1.0
    for step in range(1, 1001):
        m = step // 2
        if step % 2:
            coefficient = -(a + m) * (a + b + m) * x / ((a + 2 * m) * (a + 2 * m + 1))
        else:
            coefficient = m * (b - m) * x / ((a + 2 * m - 1) * (a + 2 * m))
        denominator = 1 + coefficient * previous_denominator
        numerator = 1 + coefficient / previous_numerator
        previous_denominator = 1 / (denominator or floor)
        previous_numerator = numerator or floor
        change = previous_denominator * previous_numerator
        denominator_fraction *= change
        if step % 2 and step > 1 and abs(change - 1) <= 1e-14:
            log_weight = a * math.log(x) + b * math.log1p(-x) - math.log(a)
            log_weight += math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
            return math.exp(log_weight) / denominator_fraction
    raise ValueError("Welch 概率计算未收敛，请核对批次电耗")


def _student_t_upper(t, degrees):
    if t == 0:
        return 0.5
    distance = abs(t)
    root_df = math.sqrt(degrees)
    norm = math.hypot(root_df, distance)
    x = (root_df / norm) ** 2
    y = (distance / norm) ** 2
    a, b = degrees / 2, 0.5
    if x == 0:
        half_tail = 0.0
    elif y == 0:
        half_tail = 0.5
    elif x < (a + 1) / (a + b + 2):
        half_tail = _regularized_beta_lower(a, b, x) / 2
    else:
        half_tail = (1 - _regularized_beta_lower(b, a, y)) / 2
    return half_tail if t > 0 else 1 - half_tail


def _upper_t_quantile(alpha, degrees):
    low, high = 0.0, 1.0
    while _student_t_upper(high, degrees) > alpha:
        high *= 2
    for _ in range(54):
        mid = (low + high) / 2
        if _student_t_upper(mid, degrees) > alpha:
            low = mid
        else:
            high = mid
    return (low + high) / 2


def analyze(rows, params):
    """Use baseline minus staged energy intensity as the positive savings direction."""
    alpha, saving, min_batches, min_inspected, max_rate, max_increase, review = _get_settings(params)
    (kiln, product), groups = _parse_batches(rows)
    control, trial = (groups[recipe] for recipe in RECIPES)
    energies = [[batch[2] for batch in group] for group in (control, trial)]
    means = [statistics.mean(sample) for sample in energies]
    deviations = [statistics.stdev(sample) if len(sample) > 1 else None for sample in energies]
    inspected = [sum(batch[3] for batch in group) for group in (control, trial)]
    defective = [sum(batch[4] for batch in group) for group in (control, trial)]
    rates = [Fraction(bad * 100, total) for bad, total in zip(defective, inspected)]
    increase = rates[1] - rates[0]
    difference = means[0] - means[1]

    coverage_ok = min(inspected) >= min_inspected
    defect_rate_ok = rates[1] <= max_rate
    increase_ok = increase <= max_increase
    p_value = t_value = degrees = lower = None
    energy_detail = f"右尾检验 H₀：平均节能量 ≤ {saving:g} kWh/千㎡，α = {alpha:g}。"
    if review == REVIEW_DEPENDENT:
        energy_reason = "存在配对或相关批次，当前数据暂停 Welch 推断。"
    elif min(len(control), len(trial)) < min_batches:
        energy_reason = f"至少一组不足 {min_batches} 个独立批次。"
    elif not any(deviations):
        energy_reason = "两组批间标准差均为零，请复核电表精度或补采批次。"
    else:
        se_control, se_trial = (sd / math.sqrt(len(sample)) for sd, sample in zip(deviations, energies))
        standard_error = math.hypot(se_control, se_trial)
        if standard_error == 0:
            raise ValueError("Welch 标准误超出可计算范围，请核对电量和面积的单位")
        # Ratios keep the Satterthwaite formula stable without squaring large variances.
        control_term = (se_control / standard_error) ** 4 / (len(control) - 1)
        trial_term = (se_trial / standard_error) ** 4 / (len(trial) - 1)
        degrees = 1 / (control_term + trial_term)
        t_value = (difference - saving) / standard_error
        lower = difference - _upper_t_quantile(alpha, degrees) * standard_error
        if not math.isfinite(lower) or not math.isfinite(t_value):
            raise ValueError("Welch 结果超出可计算范围，请核对电耗单位和最低节能量")
        p_value = _student_t_upper(t_value, degrees)
        energy_reason = ""

    reasons = []
    if not defect_rate_ok:
        reasons.append(f"分段缺陷率超过 {float(max_rate):g}%")
    if not increase_ok:
        reasons.append(f"缺陷率涨幅超过 {float(max_increase):g} 个百分点")

    if reasons:
        decision = f"维持现行预热，交质检复核（{'、'.join(reasons)}）。"
    elif review == REVIEW_DEPENDENT:
        decision = "维持现行预热，按独立批次重新安排对照后再评估。"
    elif not coverage_ok:
        decision = f"暂缓试行，补齐每组至少 {min_inspected} 片的外观抽检。"
    elif energy_reason:
        decision = f"暂缓试行。{energy_reason}"
    elif p_value > alpha:
        decision = f"维持现行预热。当前证据不足以确认节能超过 {saving:g} kWh/千㎡。"
    elif review == REVIEW_PENDING:
        decision = "暂缓试行，待烧成工艺组完成独立性、工况和分布核对。"
    else:
        decision = f"建议在 {kiln} 窑线的 {product} 规格继续限线试行分段预热，由质检与窑炉班组会签。"

    observed = f"平均节能 {difference:.2f}" if difference >= 0 else f"分段组批均电耗反升 {-difference:.2f}"
    confidence = f"{100 * (1 - alpha):g}%"
    summary = (
        f"{decision} 现行 {len(control)} 批、分段 {len(trial)} 批，"
        f"{observed} kWh/千㎡。"
        f"缺陷率由 {float(rates[0]):.3f}% 变为 {float(rates[1]):.3f}%（{float(increase):+.3f} 个百分点）。"
    )
    if p_value is not None:
        summary += f" Welch 单侧 p = {p_value:.5g}，{confidence} 单侧置信下限 {lower:.2f} kWh/千㎡。"

    def show(value, digits):
        return round(value, digits) if value is not None else "—"

    def status(passed):
        return "通过" if passed else "未通过"

    return {
        "summary": summary,
        "metrics": [
            {"label": "平均节能（kWh/千㎡）", "value": round(difference, 3), "detail": "现行批均电耗减分段批均电耗。"},
            {"label": "Welch 单侧 p 值", "value": p_value if p_value is not None else "未计算",
             "detail": energy_detail + (energy_reason or f"t = {t_value:.5g}，自由度 {degrees:.5g}。")},
            {"label": "单侧节能置信下限（kWh/千㎡）", "value": show(lower, 3),
             "detail": f"置信水平 {confidence}。"},
            {"label": "缺陷率变化（百分点）", "value": round(float(increase), 4),
             "detail": f"分段组减现行组。涨幅护栏{status(increase_ok)}。"},
        ],
        "tables": [
            {"title": f"{kiln} · {product} 对照概览",
             "columns": ["预热工艺", "录入批数", "批均电耗（kWh/千㎡）", "批间标准差", "抽检片数", "缺陷片数", "缺陷率（%）"],
             "rows": [[recipe, len(groups[recipe]), round(avg, 3), show(sd, 3), total, bad, round(float(rate), 4)]
                      for recipe, avg, sd, total, bad, rate in zip(RECIPES, means, deviations, inspected, defective, rates)]},
            {"title": "试行会签条件",
             "columns": ["条件", "规则", "实测 / 状态", "判定"],
             "rows": [
                 ["试验与分布核对", "人工确认", review, status(review == REVIEW_CONFIRMED)],
                 ["Welch 节能证据", f"节能超过 {saving:g} kWh/千㎡，p ≤ {alpha:g}",
                  p_value if p_value is not None else energy_reason,
                  status(p_value <= alpha) if p_value is not None else "待取得证据"],
                 ["抽检覆盖量（片）", f"每组 ≥ {min_inspected}", f"现行 {inspected[0]} / 分段 {inspected[1]}", status(coverage_ok)],
                 ["分段缺陷率（%）", f"≤ {float(max_rate):g}", round(float(rates[1]), 4), status(defect_rate_ok)],
                 ["缺陷率涨幅（百分点）", f"≤ {float(max_increase):g}", round(float(increase), 4), status(increase_ok)],
             ]},
            {"title": "逐批核对明细",
             "columns": ["生产日期", "批号", "预热工艺", "电耗（kWh/千㎡）", "抽检片数", "缺陷片数"],
             "rows": sorted([day, batch_id, recipe, round(intensity, 3), total, bad]
                            for recipe in RECIPES for day, batch_id, intensity, total, bad in groups[recipe])},
        ],
        "notes": [
            "电量按本批面积换算到每千平方米，组内每批等权。Welch 使用各组独立方差与 Welch–Satterthwaite 自由度。",
            "小样本推断要求批次独立且两组电耗强度近似正态。核对重复投料、窑炉余热残留、固定配方班次、明显偏态和极端批。待核对时结果仅供复核。",
            "抽检口径需一致。总缺陷片数除以总抽检片数只用于质量护栏，有限抽检无法确认两种工艺质量等效。判定使用舍入前的数值。",
            "在看数据前约定门槛和显著性水平。对照仅覆盖当前窑线和砖种，换线或换规格需重新验证。",
        ],
        "series_title": "两组批均预热电耗（kWh/千㎡）",
        "series": [{"label": f"{recipe} · 批均电耗", "value": round(avg, 3)} for recipe, avg in zip(RECIPES, means)],
    }
