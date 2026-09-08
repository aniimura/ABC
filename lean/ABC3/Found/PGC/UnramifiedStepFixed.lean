import ABC3.Found.PGC.UnramifiedLayerFree
import ABC3.Found.PGC.WildDescentDistanceOnly

/-!
# [pGC] ①の (b) が閉じた —— 不分岐の 1 段は**損失も `ε` の伸びも `1`**（完全に無料）

## ★どちらを選んだか、なぜか（持ち場が判断を任せた点）

★**(b) を選んだ。** (a) の費用を先に測った結果である:

```
grep -in "IsDiscreteValuationRing.*IsDedekindDomain|instIsDedekindDomain" .cache/mathlib-index.txt
  → ★直接の instance は無い。在るのは
    IsDiscreteValuationRing.TFAE (RingTheory/DiscreteValuationRing/TFAE.lean:210)
    —— TFAE のリストから取り出す形で、[IsNoetherianRing][IsLocalRing][IsDomain] と
       ¬IsField R が要る
grep -rn "IsDedekindDomain" lean/ABC3/Found/PGC/ --include=*.lean
  → ★木の Found/PGC には宣言が無い（UnramifiedExtension.lean:218 / :381 の docstring が
     「ここから従う」と書いているだけである）
```

⇒ ★(a) は「TFAE の取り出し」＋`Algebra.trace_quotient_eq_of_isDedekindDomain` の
残り 4 仮説（`IsDomain S` / `Module.IsTorsionFree R S` / `Module.Finite R S` /
`IsIntegrallyClosed S`）で**重い**。
★(b) は `Finset.smul_sum` と左移動の再添字だけで**軽い**と見た。★実際 1 往復で通った。

## ★★★成果 —— 不分岐の 1 段は完全に無料

§2 `exists_fixed_step_free`: `‖a‖ ≤ 1` と `Σ_{g ∈ H} g·a = 1` と
`∀ g : G, ‖g·x − x‖ ≤ ε` から

★★**`∃ y`、`H` で固定され、`‖x − y‖ ≤ ε`（損失 `1`）、かつ
`∀ g : G, ‖g·y − y‖ ≤ ε`（`ε` の伸び `1`）。**

* §1 `smul_sum_subgroup` —— `Σ_{g∈H} g·b` は `H` で固定される
  （`Finset.smul_sum` ＋ 左移動 `g ↦ h·g` の全単射）。★これが (b) の中身である。
* §1 `exists_fixed_norm_sub_le` —— 損失 `ε`
  （前波の「`‖a‖ ≤ 1` なら損失 `1`」を部分群の上の和で書き直した形）。
* `ε` の伸びは `WildDescentDistanceOnly.lean:161 norm_smul_sub_self_le_max` が
  無償で落とす（`max ε ε = ε`）。

⇒ ★★★**不分岐の層は `c = g = 1` であり、`AxTowerDecay.lean:327
exists_mem_of_descent_budget` の予算にも `AxEpsilonDecay.lean:440
axLemma_of_wildDescent_Icc` の積にも一切効かない。**

## ★残る 1 点は (a) だけ

★`‖a‖ ≤ 1` かつ `Σ_{g∈H} g·a = 1` なる `a` の**存在**（＝「不分岐 ⇒ 跡が整数環に全射」）。
★道具は mathlib に在る（前波で測った `Algebra.trace_quotient_eq_of_isDedekindDomain`、
`RingTheory/Trace/Quotient.lean:90`）が、上のとおり `IsDedekindDomain` を DVR から
出す配管が重い。★**本波では払っていない。**

## ★穴の現状

①不分岐 —— ★**(b) は閉じた**。残るのは (a) だけである。
前波の「難しい側ではない」に加えて、本波で★**「無料である」ところまで測れた**。
★ただし (a) が残るので★**「①が閉じた」とは書かない**。
②′（`ExitLossDepth`）/ ③（`JumpGeometricDecay`）/ ⑤（`ExitLossDepth`）は変わらず。④は消えたまま。

## 逸脱の記録（CLAUDE.md「逸脱」）

1. `.src` は上流と同じ pGC 物理 p.6 Corollary 3.1。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
3. 本ファイルは全部が抽象核（超距離な可換環＋等長な群作用）で、
   ★分岐・付値・Galois・`p` 進の語彙が 1 語も出ない。
4. ★`y` が「下の体に入る」ことは★**`H` で固定される**という形で書いた。
   ★`wildDepth` が真に下がることまでは示していない
   （木の `WildDepthFieldDescent.lean` の配管が要る。★本波では繋いでいない）。
-/

namespace ABC3.Found.PGC

namespace UnramifiedStep

/-! ## §1 抽象核 —— 部分群の上の和は固定され、損失は `ε` -/

section Core

variable {R : Type*} [NormedCommRing R] [IsUltrametricDist R]

omit [IsUltrametricDist R] in
/-- ★★**①の (b) の中身** —— `Σ_{g ∈ H} g·b` は `H` で固定される。

`Finset.smul_sum` で `h` を中に入れ、左移動 `g ↦ h·g` の全単射で再添字するだけ。
★分岐・付値・Galois・`p` 進の語彙が 1 語も出ない。 -/
theorem smul_sum_subgroup {G : Type*} [Group G] [MulSemiringAction G R] (H : Subgroup G)
    [Fintype H] (b : R) (h : H) :
    (h : G) • (∑ g : H, (g : G) • b) = ∑ g : H, (g : G) • b := by
  rw [Finset.smul_sum]
  refine Fintype.sum_bijective (fun g : H => h * g) (Group.mulLeft_bijective h)
    (fun g : H => (h : G) • ((g : G) • b)) (fun k : H => (k : G) • b) ?_
  intro g
  rw [← mul_smul]
  rfl

/-- ★★`‖a‖ ≤ 1`（`a` が整数）なら、`H` で固定される `y` に損失 `ε`（定数 `1`）で届く。 -/
theorem exists_fixed_norm_sub_le {G : Type*} [Group G] [MulSemiringAction G R]
    (H : Subgroup G) [Fintype H]
    (hiso : ∀ (g : G) (r : R), ‖g • r‖ = ‖r‖)
    {a x : R} {ε : ℝ} (hε : 0 ≤ ε) (ha : ‖a‖ ≤ 1)
    (hsum : ∑ g : H, (g : G) • a = 1)
    (hx : ∀ g : H, ‖x - (g : G) • x‖ ≤ ε) :
    ∃ y : R, (∀ h : H, (h : G) • y = y) ∧ ‖x - y‖ ≤ ε := by
  refine ⟨∑ g : H, (g : G) • (a * x), fun h => smul_sum_subgroup H _ h, ?_⟩
  have hrw : ∀ g : H, (g : G) • (a * x) = ((g : G) • a) * ((g : G) • x) :=
    fun g => smul_mul' (g : G) a x
  rw [Finset.sum_congr rfl (fun g _ => hrw g)]
  have h := TraceCore.norm_sub_sum_mul_le (R := R) (ι := ↥H) (s := Finset.univ)
    (a := fun g : H => (g : G) • a) (z := fun g : H => (g : G) • x) (x := x)
    (M := 1) (ε := ε) zero_le_one hε (by simpa using hsum)
    (fun g _ => le_trans (le_of_eq (hiso (g : G) a)) ha)
    (fun g _ => hx g)
  simpa using h

end Core

/-! ## §2 ★★★不分岐の 1 段は完全に無料（損失も伸びも `1`） -/

section Free

variable {R : Type*} [NormedCommRing R] [IsUltrametricDist R]

/-- ★★★★**本ファイルの主結果** —— 不分岐の 1 段は**完全に無料**。

`H` で固定される `y` が取れて、★損失 `‖x − y‖ ≤ ε` と
★`ε` の伸び `∀ g : G, ‖g·y − y‖ ≤ ε` が**両方とも定数 `1`** である。
（伸びは `WildDescentDistanceOnly.lean:161 norm_smul_sub_self_le_max` が無償で落とす。）

⇒ ★不分岐の層は予算にも積にも一切効かない。
★残るのは `a` の存在（＝不分岐 ⇒ 跡が整数環に全射）だけである。 -/
theorem exists_fixed_step_free {G : Type*} [Group G] [MulSemiringAction G R]
    (H : Subgroup G) [Fintype H]
    (hiso : ∀ (g : G) (r : R), ‖g • r‖ = ‖r‖)
    {a x : R} {ε : ℝ} (hε : 0 ≤ ε) (ha : ‖a‖ ≤ 1)
    (hsum : ∑ g : H, (g : G) • a = 1)
    (hxG : ∀ g : G, ‖g • x - x‖ ≤ ε) :
    ∃ y : R, (∀ h : H, (h : G) • y = y) ∧ ‖x - y‖ ≤ ε
      ∧ ∀ g : G, ‖g • y - y‖ ≤ ε := by
  obtain ⟨y, hfix, hd⟩ := exists_fixed_norm_sub_le H hiso hε ha hsum
    (fun g => by
      rw [norm_sub_rev]
      exact hxG (g : G))
  refine ⟨y, hfix, hd, fun g => ?_⟩
  have h := norm_smul_sub_self_le_max (G := G) (M := R) hiso hxG hd g
  simpa using h

end Free

/-! ## §3 `.src` と 使っている公理の一覧 -/

section Src

def smul_sum_subgroup.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def exists_fixed_norm_sub_le.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def exists_fixed_step_free.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end Src

#print axioms smul_sum_subgroup
#print axioms exists_fixed_norm_sub_le
#print axioms exists_fixed_step_free

end UnramifiedStep

end ABC3.Found.PGC
