import ABC3.Found.PGC.LocalFieldNorm

/-!
# Track B — log-shell は構成できるか(2026-08-14 の実測)

## 原文の定義(目視確認済み)

原典: S. Mochizuki, *Topics in Absolute Anabelian Geometry III* [AbsTopIII]、
物理 p.3(`0_Source/Topics in Absolute Anabelian Geometry III.pdf`、
全 164 ページ。**400 dpi 目視確認 2026-08-14**)。

原文 (AbsTopIII p.3):
> One important property of the p-adic logarithm log

原文 (AbsTopIII p.3):
> — which is compact — may be thought of as defining a sort of canonical rigid
> integral structure on k. In the present paper, we shall refer to the “canonical
> rigid integral structures” obtained in this way as log-shells.

その間に置かれた表示式は(400 dpi 目視):

```
log_{k[bar]}(O[scr]^×_k)  ⊆  k
```

すなわち **log-shell とは「整数環の単数群 `𝒪_k^×` の、p進対数による像」**であり、
それが**コンパクト**である、というのが原文の主張。
★`log` の添字 `k̄` にはオーバーバーが付く(pdftotext は落とす。IUTchIII で
何度も踏んだのと同じ罠。ここでは `[bar]` で明示した)。

なお [AbsTopIII] 物理 p.5 の (L1) も同じことを述べる:

原文 (AbsTopIII p.5):
> (L1) a log-shell is compact and hence of finite “log-volume” [cf. Corollary

## ★★結論 — **log-shell は構成できなかった。理由は1点に絞れる。**

**p進対数が mathlib に無い。** 探索範囲は下記(`LeanStatus.absent (searched)` の規律、
`tools/check.mjs` 冒頭 A2 の手順 S1–S4)。

- **(S1) 原文側の呼称**: 「the **p-adic logarithm** `log_{k̄}`」——
  単数群 `𝒪_k^×`(あるいは `k^×`)から加法群 `k` への写像。
- **(S2) 全出現の列挙**:
  - `Mathlib/NumberTheory/Padics/` は全 13 項目
    (AddChar, Complex, HeightOneSpectrum, Hensel, MahlerBasis, PadicIntegers,
    PadicNorm, PadicNumbers, PadicVal, ProperSpace, RingHoms, ValuativeRel, WithVal)。
    このうち名前に `log` を含む宣言は 5 件で、**すべて `Nat.log` / 付値の対数**
    (`norm_eq_zpow_log_mulValuation`, `padicValNat_le_nat_log` 等)。
  - `grep -rniE "p-adic (logarithm|exponential)" Mathlib/` → **0 件**。
  - `Mathlib/Analysis/SpecialFunctions/Log/` は 12 ファイル。すべて
    実・複素・`ENNReal`・`EReal` 用(`Real.log`, `Real.logb`)。
    **ノルム環 / Banach 環用の `log` は無い**(`NormedSpace.log` は存在しない)。
  - `NormedSpace.exp`(`Analysis/SpecialFunctions/Exponential.lean`)は**あるが**、
    対応する `log` は無い。
- **(S3) Lean 側の語で絞る**: `padicLog` / `Padic.log` / `padic_log` /
  `logOneAdd` / `log_one_add` → 6 件、すべて実/複素の解析。
  `padicExp` / `expPadic` → 0 件。
- **(S4) 記録**: 上記のとおり。**この測定は
  `ResearchPaper/lean-ecosystem.json` の `pGC#6 p進対数` = `absent` と同じ穴である**
  ——pGC Proposition 1.2 と IUT の log-shell を**同一の欠落**が塞いでいる。

## 到達点 — ここまでは `sorry` 無しで作れた

log-shell そのものは作れないが、**その定義域**は作れる。しかも
コンパクト性の証明は「コンパクト集合の連続像はコンパクト」だけなので、
**p進対数さえあれば log-shell のコンパクト性は1行で出る**。それを下で確定させる:

1. `coe_units_integer_eq_sphere` — `𝒪[K]` の**単数群**は、K の中では
   ちょうど**単位球面** `‖x‖ = 1` である。
2. `isCompact_unitSphere` — その球面は**コンパクト**(`ProperSpace K` から)。
3. `isCompact_logShell_of_continuous` — ★**連続写像 `f` を与えれば、
   `f '' (単数球面)` はコンパクト**。`f` に p進対数を入れれば、これが log-shell の
   コンパクト性そのものになる。

**3 は log-shell の構成ではない。** `f` は仮引数であって、我々は p進対数を作っていない。
ここを取り違えないこと——`f := id` と置けば「log-shell」は単数球面になってしまう
(退化)。退化と本物を分けるのは `f` が**本当に p進対数であること**だけで、
それが今まさに作れないものである。
-/

namespace ABC3.Found.IUTchIII

open scoped NormedField Valued

variable {K : Type*} [NontriviallyNormedField K] [IsUltrametricDist K]

set_option linter.unnecessarySimpa false in
/-- `𝒪[K]` の単元は、`K` の中ではちょうど**ノルム 1 の元**。

`⊆` は `Valuation.Integers.one_of_isUnit`、`⊇` は
`Valuation.Integer.not_isUnit_iff_valuation_lt_one` による。 -/
theorem coe_units_integer_eq_sphere :
    (fun u : (𝒪[K])ˣ => ((u : 𝒪[K]) : K)) '' Set.univ = Metric.sphere (0 : K) 1 := by
  ext x
  simp only [Set.image_univ, Set.mem_range, Metric.mem_sphere, dist_zero_right]
  constructor
  · rintro ⟨u, rfl⟩
    have := (Valuation.integer.integers (Valued.v (R := K))).one_of_isUnit u.isUnit
    simpa [← NNReal.coe_inj] using this
  · intro hx
    have hmem : x ∈ 𝒪[K] := by
      rw [Valued.integer.mem_iff]
      simpa [← NNReal.coe_le_coe, hx] using le_refl (1 : ℝ)
    refine ⟨(IsUnit.unit (?_ : IsUnit (⟨x, hmem⟩ : 𝒪[K]))), ?_⟩
    · by_contra hnu
      rw [Valuation.Integer.not_isUnit_iff_valuation_lt_one] at hnu
      have : ‖x‖ < 1 := by simpa [← NNReal.coe_lt_coe] using hnu
      exact absurd hx (ne_of_lt this)
    · rfl

variable [ProperSpace K]

/-- 単数球面は**コンパクト** —— `K` が proper だから。 -/
theorem isCompact_unitSphere : IsCompact (Metric.sphere (0 : K) 1) :=
  isCompact_sphere 0 1

/-- ★**p進対数さえあれば、log-shell のコンパクト性はこれだけ**。

原文 (AbsTopIII p.5):
> (L1) a log-shell is compact and hence of finite “log-volume” [cf. Corollary

★これは **log-shell の構成ではない**。`f` は仮引数であり、我々は p進対数を
作っていない(`f := id` と置けば退化する)。この定理が言うのは
「欠けているのは `f` だけで、コンパクト性の議論そのものは既に手元にある」ということ。 -/
theorem isCompact_logShell_of_continuous {f : K → K} (hf : Continuous f) :
    IsCompact (f '' Metric.sphere (0 : K) 1) :=
  (isCompact_sphere 0 1).image hf

/-! ## p進局所体への特殊化 -/

open ABC3.Skeleton.PGC ABC3.Found.PGC

variable {p : ℕ} [Fact p.Prime]

/-- `PAdicLocalField p` の単数球面はコンパクト。 -/
theorem isCompact_unitSphere_padic (K : PAdicLocalField p) :
    IsCompact (Metric.sphere (0 : K.carrier) 1) :=
  isCompact_sphere 0 1

/-- `PAdicLocalField p` について、p進対数(連続写像)があれば log-shell はコンパクト。 -/
theorem isCompact_logShell_padic (K : PAdicLocalField p)
    {f : K.carrier → K.carrier} (hf : Continuous f) :
    IsCompact (f '' Metric.sphere (0 : K.carrier) 1) :=
  (isCompact_sphere 0 1).image hf

end ABC3.Found.IUTchIII
