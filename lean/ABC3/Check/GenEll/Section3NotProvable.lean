import ABC3.Interface.GenEll.TateLocal
import ABC3.Interface.GenEll.EllModuli
import ABC3.Found.GenEll.BDClass

/-!
# 退化封じの検査 —— **§3 の 5 つの `sorry` は「難しい」のではなく偽である**(`Check`)

**これは原典の主張ではない**(我々のモデルについての事実)ので `.src` を持たない。

## ★★★★★★★★2026-08-26 の診断

`Skeleton/GenEll/Section3.lean` に残る 5 つの `sorry` は

| 定理 | 量化している界面 |
|---|---|
| `lemma_3_2` | `∀ D : TateLocalData` |
| `potLocalHeight_indep` | `∀ D : TateLocalData` ★**2026-08-26 に閉じた** |
| `prop_3_4` | `∀ D : EllModuliData` |
| `lemma_3_5` | `∀ D : EllModuliData` |
| `lemma_3_7` | `∀ D : EllModuliData` |

という形をしている。★★**ところがこの 2 つの界面は公理を 1 つも持たない**
(`TateLocalData` は `vq_pos`・`ramIdx_pos` の 2 本、`EllModuliData` は 0 本)。
★★★したがって**退化した `D` を作れば主張は破れる**——`sorry` は
「まだ証明していない」ではなく「**この形では証明できない**」ことの印である。

## ★★★★★★これは G6 が 2026-08-17 に通った道と同じである

`Interface/GaloisRep/Reduction.lean` の `TateCurveData` は、以前
`LocalField : Type` / `Curve : LocalField → Type` と**世界ごと posit** していて
`PUnit` と定数 `1` で埋まった。★そこで mathlib の `WeierstrassCurve` と
正規化付値に**接地**し、**Tate 一意化そのもの**を要求する形に作り直した。

★★`Interface/GenEll/TateLocal.lean` は**その捨てたはずの形のまま**である
(`LocalField : Type` / `Curve : LocalField → Type`)。
★★★`EllModuliData` も `EllClass : Type` / `faltingsHeight : EllClass → ℝ` を
条件なしで持つだけである。

## ★★★★本ファイルが示すこと(限界の明示)

★示すのは「**現在の界面の上では §3 の 5 つは偽である**」ことだけである。
★★「原文の `Lemma 3.2` / `Proposition 3.4` が偽」ということでは**まったくない**。
★★★塞ぎ方も分かっている——G6 と同じく、**既に達成済みの
`TateCurveData`(G6)・`FaltingsHeightData`(G8)に接地する**ことである。

## ★★界面の欠陥の一覧(本ファイルで 2 つ増える)

| # | 場所 | 欠陥の型 | 塞いだ |
|---|---|---|---|
| 1-5 | G6/G7/G8 | 充足不能・弱すぎる | 第 302-318 |
| 6 | G8 `htFalt` | **弱すぎる**(`deg∞/12` で満たせる) | 第 357 |
| 7 | **`TateLocalData`** | **弱すぎる**(世界ごと posit、`Unit` で埋まる) | ★**部分**(第 360: `Remark 3.3.1`) |
| 8 | **`EllModuliData`** | **弱すぎる**(高さの間に関係が無い) | ★**未** |
-/

namespace ABC3.Check.GenEll

open ABC3.Interface.GenEll ABC3.Found.GenEll

/-! ## ★★★★★★1. `TateLocalData` は `Unit` で埋まる -/

/-- ★**退化した `TateLocalData`**——世界を `Unit`、`v_K(q_E)` を定数 `1` にしたもの。

★`MultExt` だけは 2 元にしてある(`Remark 3.3.1` を破るため)。 -/
def degenerateTateLocal : TateLocalData where
  LocalField := Unit
  Curve := fun _ => Unit
  vq := fun _ _ => 1
  vq_pos := fun _ _ => Nat.one_pos
  degInf := fun _ _ => 1
  StableLine := fun _ _ _ => Unit
  IsCyclotomic := fun _ => False
  quotMu := fun _ _ _ => ()
  MultExt := fun _ _ => Bool
  extField := fun _ => ()
  baseChange := fun _ => ()
  ramIdx := fun _ => 1
  ramIdx_pos := by
    intro K E L
    norm_num
  vq_baseChange := by
    intro K E L L'
    rfl

/-- ★★★★★★**`Lemma 3.2, (i)` は現在の界面の上では偽である**。

`l = 2`、`v_K(q_E) = 1`、`IsCyclotomic := False` で破れる。 -/
theorem lemma_3_2_i_false :
    ¬ (∀ (K : degenerateTateLocal.LocalField) (E : degenerateTateLocal.Curve K) (l : ℕ)
        (N : degenerateTateLocal.StableLine K E l),
      l ∣ degenerateTateLocal.vq K E ∨ degenerateTateLocal.IsCyclotomic N) := by
  intro h
  rcases h () () 2 () with hdvd | hcyc
  · exact absurd hdvd (by decide)
  · exact hcyc

/-- ★★★★★★**`Lemma 3.2, (ii)` は現在の界面の上では偽である**。

`deg_∞` が定数なので `deg_∞(E/μ_l) = l·deg_∞(E)` は `l = 2` で破れる。 -/
theorem lemma_3_2_ii_false :
    ¬ (∀ (K : degenerateTateLocal.LocalField) (E : degenerateTateLocal.Curve K) (l : ℕ),
      degenerateTateLocal.degInf K (degenerateTateLocal.quotMu K E l)
        = (l : ℝ) * degenerateTateLocal.degInf K E) := by
  intro h
  have h2 := h () () 2
  show False
  simp only [degenerateTateLocal] at h2
  norm_num at h2

/-! ## ★★★★★★★`Remark 3.3.1` は 2026-08-26 に閉じた

★以前はここに `potLocalHeight_indep_false` があった——分岐指数だけが違う
2 つの拡大を取れば `v/e` は一致せず、`Remark 3.3.1` は破れていた。

★★`Interface/GenEll/TateLocal.lean` に **`vq_baseChange`**(局所高さは分岐指数倍)を
足したので、**その退化 witness はもはや作れない**。
★★★これは posit ではなく定理である——`Found/GenEll/LocalHeightRamified.lean`(第 359)が
mathlib の `HeightOneSpectrum.valuation_liesOver` から導いている。
★★★★`Skeleton/GenEll/Section3.lean` の `potLocalHeight_indep` はこれで
**`sorry` ではなくなった**。

★下の `degenerateTateLocal` は `ramIdx ≡ 1` にしてあるので新しい欄を満たすが、
`Lemma 3.2` は依然として破る。 -/

/-! ## ★★★★★★2. `EllModuliData` は高さの間に関係を持たない -/

/-- ★**退化した `EllModuliData`**——`ht_∞` だけを非有界にし、他をすべて自明にしたもの。 -/
def degenerateEllModuli : EllModuliData where
  EllClass := ℕ
  Curve := Unit
  cls := fun _ => 0
  degOfDefinition := fun _ => 1
  faltingsHeight := fun _ => 0
  HasPotMultRed := fun _ => True
  PrimeToLocalHeights := fun _ _ => True
  CompactlyBounded := fun _ => True
  GaloisFinite := fun _ => True
  ImageContainsSL2 := fun _ _ => True
  degInf := fun _ => 0
  htInf := fun n => (n : ℝ)
  logDiffMell := fun _ => 0
  degLe := fun _ => Set.univ
  SemiStable := fun _ => True
  HasLCyclic := fun _ _ => True
  MinimalField := fun _ => True
  ImageSurjective := fun _ _ => True
  PrimeToRamification := fun _ _ => True
  HasMultRed := fun _ => True
  PrimeToMultPrimes := fun _ _ => True

/-- ★★★★★★**`Proposition 3.4` の最初の `≲` は現在の界面の上では偽である**。

`deg_∞ ≡ 0`・`ht_∞ = n` とすれば `ht_∞ − deg_∞` は非有界。 -/
theorem prop_3_4_first_false :
    ¬ BDle degenerateEllModuli.degInf degenerateEllModuli.htInf := by
  rintro ⟨C, hC⟩
  obtain ⟨n, hn⟩ := exists_nat_gt C
  have h2 := hC n
  simp only [degenerateEllModuli, sub_zero] at h2
  linarith

/-- ★★★★★★**`Proposition 3.4` の有限性(Northcott)も現在の界面の上では偽である**。

`ht^Falt ≡ 0`・`M_ell^{≤d} = univ` なので、集合は `ℕ` 全体になる。 -/
theorem prop_3_4_finite_false :
    ¬ (∀ (C : ℝ) (d : ℕ), 0 < d →
      {x ∈ degenerateEllModuli.degLe d | degenerateEllModuli.faltingsHeight x ≤ C}.Finite) := by
  intro h
  have h2 := h 0 1 Nat.one_pos
  haveI : Infinite degenerateEllModuli.EllClass := (inferInstance : Infinite ℕ)
  have hsub : (Set.univ : Set degenerateEllModuli.EllClass)
      ⊆ {x ∈ degenerateEllModuli.degLe 1 | degenerateEllModuli.faltingsHeight x ≤ 0} := by
    intro n _
    exact ⟨trivial, le_refl 0⟩
  exact (Set.Infinite.mono hsub Set.infinite_univ) h2

/-! ## ★★★塞ぎ方は分かっている -/

/-- ★★★★**`EllModuliData` は `ht^Falt ≡ 0` を許す**——G8 で本物の `ht^Falt` を
建てた(第 357)ので、そこに接地すれば塞がる。 -/
theorem faltingsHeight_can_be_zero :
    ∃ D : EllModuliData, ∀ x : D.EllClass, D.faltingsHeight x = 0 :=
  ⟨degenerateEllModuli, fun _ => rfl⟩

end ABC3.Check.GenEll
