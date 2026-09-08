import ABC3.Found.PGC.NormalizedTraceDescent

/-!
# [pGC] ★①不分岐側は「難しい側」ではなく**最良定数 `1` の側**である

## ★どちらを選んだか、なぜか（持ち場が判断を任せた点）

★**①（不分岐側）を選んだ。** 理由は 3 つ:

1. 前波で★**初めて①が本筋に入った**（`FirstJumpRoute` の `heq : e_{E₁} = p·e_K` は
   完全分岐の帰結で、不分岐なら偽）。持ち場も「先に測る価値がある」と書いた。
2. `hform`（第 1 跳びの損失を達成する降下の構成）は★**新しい降下の機構**が要ると測れた
   （前波の報告のとおり、`lt_ramIndex_quotient_fixedRing` は「`σ̄` が `x` を動かす量」の
   下界であって「`x′` の作り方」ではない）。★1 波では入らないと見た。
3. ★①は**本日ずっと未着手**で、しかも 4 か所の docstring が「別の議論が要る」と
   書いているだけで、**中身を測った波が 1 つも無い**。

## ★★★測定の結論 —— ①は「別の議論」ではなく「**無料の側**」である

跡による降下（`NormalizedTraceDescent.lean:179 TraceCore.norm_sub_traceAverage_le`）の
損失は `‖a‖`（`a` は `Σ_{g∈s} g·a = 1` を満たす元）である。本ファイルは:

* §1 `one_le_norm_of_sum_smul_eq_one` —— ★★**`1 ≤ ‖a‖` は無条件**
  （超距離で `1 = ‖Σ g·a‖ ≤ max_g ‖g·a‖ = ‖a‖`）。
  ⇒ ★**跡による降下の損失は決して `1` を下回れない。`1` が最良定数である。**
* §1 `norm_sub_traceAverage_le_of_norm_le_one` —— ★`‖a‖ ≤ 1` なら損失は **`ε`（定数 `1`）**。
* §1 `norm_eq_one_of_sum_smul_eq_one` —— したがって `‖a‖ ≤ 1` は `‖a‖ = 1` と同値。
* §2 `traceAverage_best_constant` —— 上の 2 つを 1 本にした形。

★**`‖a‖ ≤ 1`（＝ `a` が整数）が取れるのは「跡が整数環の上に全射」のときであり、
それはちょうど `L/M` が不分岐（different が自明）のときである。**
⇒ ★★★**不分岐の層は損失 `1` で通り、積の勘定に一切効かない。**
★分岐した層でだけ `‖a‖ > 1` を払う（`NormalizedTraceDescent.lean:262
traceLoss_eq_sharpLoss_rpow` が測った `(p−1)` 倍がその中身）。

## ★★木の断定の読み方を訂正する（名指しで書く。他ファイルは書き換えない）

`TotallyRamifiedLayer.lean:328-341` は
「★不分岐の側では `π` が `‖K^×‖` に入ってしまい、出口の `hvalK` は偽になる」
「⇒ ★不分岐の側は**別の議論**が要る」と書いている。

★**前半（`hvalK` が偽）は真である**（`not_valK_of_norm_eq`。本日 3 度読んで確かめた）。
★★しかし後半の「別の議論が要る」は★**「難しい」という意味ではない**。
上の測定のとおり、★不分岐の層は**跡で損失 `1`**、すなわち★**最良定数**で通る。
偽になるのは「全分岐を仮定した**その出口**の仮説」であって、★降下そのものではない。

★同 `:37-41` の「`p ∤ deg minpoly` と『層が不分岐』を同一視しないこと」という警告は
★**今も正しい**（本ファイルは同一視していない —— §1 は `‖a‖ ≤ 1` だけを仮定し、
次数にも分岐にも触れていない）。

## ★★まだ測っていないこと（正直に）

* ★**「不分岐 ⇒ `∃ a ∈ 𝒪_L, Tr(a) = 1`」を形式化していない。**
  本ファイルは `‖a‖ ≤ 1` を**仮説として受けている**。
  ★これは `Algebra.intTrace` の全射性（different が自明）である。
  ★★**道具は mathlib に在る（本波で 2 か所測った）**:
  ```
  grep -in "intTrace|trace_surjective" .cache/mathlib-index.txt
    → Algebra.intTrace (RingTheory/IntegralClosure/IntegralRestrict.lean:260)
    → ★Algebra.trace_quotient_eq_of_isDedekindDomain (RingTheory/Trace/Quotient.lean:90)
        Tr_{k_L/k_M}(x mod) = intTrace(x) mod
  grep -rn "intTrace" lean/ABC3/Found/PGC/ --include=*.lean
    → 木の Found/PGC には無い（GenEll/DifferentTameGlobal.lean に在ると
       NormalizedTraceDescent.lean:121 が測っている）
  ```
  ⇒ ★不分岐なら剰余体の拡大が分離的なので `Tr_{k_L/k_M}` が全射、
  よって `intTrace(x)` が単元になる `x` が取れ、`a := x / intTrace(x)` が欲しいものである。
  ★★**これが次の 1 点の正確な形である**（本波では配管が重いので形式化していない）。
* ★平均化した `Σ_{g∈s} g·(a·x)` が下の体に入ること（深さが下がること）も
  本ファイルでは示していない（木の `WildDepthFieldDescent` の配管が持っている）。

## ★穴の現状

①不分岐 —— ★★**「難しい側」ではないと測れた**（損失 `1`、最良）。
★ただし上の未測定 2 点が残るので★**閉じたとは書かない**。
②′ / ③幾何減衰 / ⑤出口の限界 は変わらず。④は消えたまま。

## 逸脱の記録（CLAUDE.md「逸脱」）

1. `.src` は `NormalizedTraceDescent.lean` と同じ項目（pGC 物理 p.6 Corollary 3.1）。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
3. 本ファイルは全部が抽象核である —— 超距離な可換環と等長な群作用だけで、
   ★分岐・付値・Galois・`p` 進の語彙が 1 語も出ない。
4. ★`‖a‖ ≤ 1` を仮説のまま残したことを上に明記した。
-/

namespace ABC3.Found.PGC

namespace UnramifiedFree

/-! ## §1 抽象核 —— 跡の損失は必ず `1` 以上、`1` は「`a` が整数」のとき -/

section Core

variable {R : Type*} [NormedCommRing R] [IsUltrametricDist R]

/-- ★★★**跡の損失は決して `1` を下回れない** —— `Σ_{g∈s} g·a = 1` なら `1 ≤ ‖a‖`。

超距離で `1 = ‖Σ g·a‖ ≤ max_g ‖g·a‖ = ‖a‖`（等長性）。
⇒ ★`1` が跡による降下の**最良定数**である。
★分岐・付値・Galois・`p` 進の語彙が 1 語も出ない。 -/
theorem one_le_norm_of_sum_smul_eq_one {G : Type*} [Group G] [MulSemiringAction G R]
    (hiso : ∀ (g : G) (r : R), ‖g • r‖ = ‖r‖) (h1 : ‖(1 : R)‖ = 1)
    {s : Finset G} {a : R} (hsum : ∑ g ∈ s, g • a = 1) : 1 ≤ ‖a‖ := by
  have hle : ‖∑ g ∈ s, g • a‖ ≤ ‖a‖ :=
    IsUltrametricDist.norm_sum_le_of_forall_le_of_nonneg (norm_nonneg a)
      (fun g _ => le_of_eq (hiso g a))
  rw [hsum, h1] at hle
  exact hle

/-- ★★**`a` が整数（`‖a‖ ≤ 1`）なら損失は `ε`（定数 `1`）**。

★`‖a‖ ≤ 1` が取れるのは「跡が整数環の上に全射」のとき、すなわち
★**層が不分岐**（different が自明）のときである。
⇒ ★不分岐の層は**最良定数 `1`** で通り、積の勘定に効かない。 -/
theorem norm_sub_traceAverage_le_of_norm_le_one {G : Type*} [Group G] [MulSemiringAction G R]
    (hiso : ∀ (g : G) (r : R), ‖g • r‖ = ‖r‖)
    {s : Finset G} {a x : R} {ε : ℝ} (hε : 0 ≤ ε) (ha : ‖a‖ ≤ 1)
    (hsum : ∑ g ∈ s, g • a = 1) (hx : ∀ g ∈ s, ‖x - g • x‖ ≤ ε) :
    ‖x - ∑ g ∈ s, g • (a * x)‖ ≤ ε := by
  have h := TraceCore.norm_sub_traceAverage_le hiso hε hsum hx
  calc ‖x - ∑ g ∈ s, g • (a * x)‖ ≤ ‖a‖ * ε := h
    _ ≤ 1 * ε := mul_le_mul_of_nonneg_right ha hε
    _ = ε := one_mul ε

/-- したがって `‖a‖ ≤ 1` は `‖a‖ = 1` と同値である。 -/
theorem norm_eq_one_of_sum_smul_eq_one {G : Type*} [Group G] [MulSemiringAction G R]
    (hiso : ∀ (g : G) (r : R), ‖g • r‖ = ‖r‖) (h1 : ‖(1 : R)‖ = 1)
    {s : Finset G} {a : R} (hsum : ∑ g ∈ s, g • a = 1) (ha : ‖a‖ ≤ 1) : ‖a‖ = 1 :=
  le_antisymm ha (one_le_norm_of_sum_smul_eq_one hiso h1 hsum)

end Core

/-! ## §2 まとめと重心版 -/

section Best

variable {R : Type*} [NormedCommRing R] [IsUltrametricDist R]

/-- ★§1 の 2 本を 1 つにした形。 -/
theorem traceAverage_best_constant {G : Type*} [Group G] [MulSemiringAction G R]
    (hiso : ∀ (g : G) (r : R), ‖g • r‖ = ‖r‖) (h1 : ‖(1 : R)‖ = 1)
    {s : Finset G} {a : R} (hsum : ∑ g ∈ s, g • a = 1) :
    1 ≤ ‖a‖ ∧ (‖a‖ ≤ 1 ↔ ‖a‖ = 1) := by
  have hge := one_le_norm_of_sum_smul_eq_one hiso h1 hsum
  exact ⟨hge, ⟨fun h => le_antisymm h hge, fun h => le_of_eq h⟩⟩

/-- 重心版（`(|s|)·a = 1`）。★`p ∤ |s|`（tame）なら `‖a‖ = 1` になる。 -/
theorem barycenter_loss_le_of_norm_le_one {ι : Type*} {s : Finset ι} {z : ι → R}
    {x a : R} {ε : ℝ} (hε : 0 ≤ ε) (hna : (s.card : R) * a = 1) (ha : ‖a‖ ≤ 1)
    (hz : ∀ i ∈ s, ‖x - z i‖ ≤ ε) :
    ‖x - ∑ i ∈ s, a * z i‖ ≤ ε := by
  have h := TraceCore.norm_sub_barycenter_le hε hna hz
  calc ‖x - ∑ i ∈ s, a * z i‖ ≤ ‖a‖ * ε := h
    _ ≤ 1 * ε := mul_le_mul_of_nonneg_right ha hε
    _ = ε := one_mul ε

end Best

/-! ## §3 `.src` と 使っている公理の一覧 -/

section Src

def one_le_norm_of_sum_smul_eq_one.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def norm_sub_traceAverage_le_of_norm_le_one.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def traceAverage_best_constant.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end Src

#print axioms one_le_norm_of_sum_smul_eq_one
#print axioms norm_sub_traceAverage_le_of_norm_le_one
#print axioms norm_eq_one_of_sum_smul_eq_one
#print axioms traceAverage_best_constant
#print axioms barycenter_loss_le_of_norm_le_one

end UnramifiedFree

end ABC3.Found.PGC
