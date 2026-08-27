import ABC3.Found.NCBelyi.ReducedDegree

/-!
# [NCBelyi] Lemma 2.4 —— 入れ子帰納法の測度(`Found`)

原典: S. Mochizuki, *Noncritical Belyi Maps* [NCBelyi]、物理 p.5。

原文 (NCBelyi p.5):
> Now we perform a "nested induction" on m(S), d(S): That is to say, we induct

## ★★★本ファイルが取るもの

原文は `m(S)`(reduced degree の最大値)と `d(S)`(最大値を取る点の reduced degree の和)
を定め、辞書式に下がることで帰納法を回す。★その **1 段が下がること**をここで取る。

原文 (NCBelyi p.5):
> m(S ) = m(S);

★★原文は結論(`m(S′) < m(S)` または `m(S′) = m(S), d(S′) < d(S)`)だけを書く。
その根拠は 3 つで、実装では**仮定として明示した**:

| 根拠 | 仮定 |
|---|---|
| 像は次数を上げない | `hno : ∀ x ∈ S, redDeg (g x) ≤ redDeg x` |
| 余分な点(`f₀(S₀)`)は最大層より下 | `hE : ∀ y ∈ E, redDeg y < m(S)` |
| 最大層の点が 1 つ落ちる(`f₀(α₀) = 0`) | `hx₀drop : redDeg (g x₀) < m(S)` |

★★★**3 つめが本質である。** 1 つめと 2 つめだけでは `m` も `d` も下がらない
——`g` が最大層を最大層へ全単射で写せば測度は不変になる。
原文が `f₀` を **`α₀` の最小多項式**に取るのは、まさに `f₀(α₀) = 0` を得るためである。

## ★★★★`d(S)` は「層の要素数 × 層の高さ」

原文の `d(S)` は「`m(S)` に等しい reduced degree の総和」だが、
★和の各項が同じ値 `m(S)` なので、**`|top(S)| · m(S)`** である(`dSum_eq`)。
★★だから `m(S) > 0` の下では `d` が下がることと**層の要素数が減ること**が同値になる。
-/

namespace ABC3.Found.NCBelyi

open Polynomial

/-! ## ★`m(S)`、`top(S)`、`d(S)` -/

/-- **`m(S)`** —— `S` の点の reduced degree の最大値。

原文 (NCBelyi p.5):
> for the maximum of the reduced degrees of the fields of definition of the various -/
noncomputable def maxRedDeg (S : Finset ℂ) : ℕ := S.sup redDeg

/-- **`top(S)`** —— reduced degree が `m(S)` に等しい点の集まり(原文の `d(S)` の台)。 -/
noncomputable def topLayer (S : Finset ℂ) : Finset ℂ :=
  S.filter (fun x => redDeg x = maxRedDeg S)

/-- **`d(S)`** —— `top(S)` の reduced degree の和。

原文 (NCBelyi p.5):
> for the sum of those reduced degrees of the fields of definition of the various points -/
noncomputable def dSum (S : Finset ℂ) : ℕ := ∑ x ∈ topLayer S, redDeg x

theorem redDeg_le_maxRedDeg {S : Finset ℂ} {x : ℂ} (hx : x ∈ S) :
    redDeg x ≤ maxRedDeg S := Finset.le_sup hx

theorem mem_topLayer_iff {S : Finset ℂ} {x : ℂ} :
    x ∈ topLayer S ↔ x ∈ S ∧ redDeg x = maxRedDeg S := Finset.mem_filter

/-- ★★**`d(S) = |top(S)| · m(S)`** —— 和の各項が同じ値だから。 -/
theorem dSum_eq (S : Finset ℂ) : dSum S = (topLayer S).card * maxRedDeg S := by
  rw [dSum, Finset.sum_congr rfl (fun x hx => (mem_topLayer_iff.1 hx).2), Finset.sum_const,
    smul_eq_mul]

/-- ★**`m(S) = 0` は「`S` がすべて有理点」と同値**。

原文 (NCBelyi p.5):
> f and only if m(S) = 0 if and only if S ⊆P1(Q). -/
theorem maxRedDeg_eq_zero_iff {S : Finset ℂ} :
    maxRedDeg S = 0 ↔ ∀ x ∈ S, redDeg x = 0 := by
  constructor
  · intro h x hx
    have := redDeg_le_maxRedDeg hx
    omega
  · intro h
    refine Nat.le_zero.1 (Finset.sup_le (fun x hx => ?_))
    exact le_of_eq (h x hx)

/-- ★★**`d(S) = 0 ⟺ m(S) = 0`**。

原文 (NCBelyi p.5):
> f and only if m(S) = 0 if and only if S ⊆P1(Q).

★`S` が空でないときの主張である(空なら `top(S)` も空で両辺 `0`)。 -/
theorem dSum_eq_zero_iff_of_nonempty {S : Finset ℂ} (hne : S.Nonempty) :
    dSum S = 0 ↔ maxRedDeg S = 0 := by
  obtain ⟨z, hzS, hzv⟩ := Finset.exists_mem_eq_sup S hne redDeg
  have hztop : z ∈ topLayer S := mem_topLayer_iff.2 ⟨hzS, hzv.symm⟩
  rw [dSum_eq]
  constructor
  · intro h
    rcases Nat.mul_eq_zero.1 h with hc | hc
    · exact absurd hc (Finset.card_ne_zero_of_mem hztop)
    · exact hc
  · intro h
    rw [h, mul_zero]

/-! ## ★★★★★★★測度が 1 段下がること -/

/-- ★★★★★★★★**入れ子帰納法の 1 段**。

原文 (NCBelyi p.5):
> m(S ) = m(S);

`T ≝ g(S) ∪ E` について、
(1) `g` が次数を上げず、(2) `E` が最大層より下にあり、
(3) 最大層の点 `x₀` が 1 つ落ちるなら、
**`m(T) < m(S)`** または **`m(T) = m(S)` かつ `d(T) < d(S)`** である。

★★★等式の場合の核心は `top(T) ⊆ g(top(S) ∖ {x₀})` である
——`redDeg (g x) = m(S)` なら `redDeg x = m(S)` が挟み撃ちで出て、
`x = x₀` は `hx₀drop` で排除される。 -/
theorem measure_lt (S E : Finset ℂ) (g : ℂ → ℂ)
    (hno : ∀ x ∈ S, redDeg (g x) ≤ redDeg x)
    (hE : ∀ y ∈ E, redDeg y < maxRedDeg S)
    {x₀ : ℂ} (hx₀S : x₀ ∈ S) (hx₀top : redDeg x₀ = maxRedDeg S)
    (hx₀drop : redDeg (g x₀) < maxRedDeg S)
    (hm0 : 0 < maxRedDeg S) :
    maxRedDeg (S.image g ∪ E) < maxRedDeg S
      ∨ (maxRedDeg (S.image g ∪ E) = maxRedDeg S ∧ dSum (S.image g ∪ E) < dSum S) := by
  classical
  set T : Finset ℂ := S.image g ∪ E with hTdef
  -- `m(T) ≤ m(S)`
  have hmle : maxRedDeg T ≤ maxRedDeg S := by
    refine Finset.sup_le (fun y hy => ?_)
    rw [hTdef, Finset.mem_union] at hy
    rcases hy with hy | hy
    · obtain ⟨x, hxS, rfl⟩ := Finset.mem_image.1 hy
      exact le_trans (hno x hxS) (redDeg_le_maxRedDeg hxS)
    · exact (hE y hy).le
  rcases lt_or_eq_of_le hmle with hlt | heq
  · exact Or.inl hlt
  refine Or.inr ⟨heq, ?_⟩
  -- `top(T) ⊆ g(top(S) ∖ {x₀})`
  have hsub : topLayer T ⊆ ((topLayer S).erase x₀).image g := by
    intro y hy
    obtain ⟨hyT, hyv⟩ := mem_topLayer_iff.1 hy
    rw [heq] at hyv
    have hyE : y ∉ E := fun hc => absurd hyv (by have := hE y hc; omega)
    rw [hTdef, Finset.mem_union] at hyT
    obtain ⟨x, hxS, rfl⟩ := Finset.mem_image.1 (hyT.resolve_right hyE)
    have hxtop : redDeg x = maxRedDeg S := by
      have h1 := hno x hxS
      have h2 := redDeg_le_maxRedDeg hxS
      omega
    have hxne : x ≠ x₀ := by
      intro hc
      rw [hc] at hyv
      omega
    exact Finset.mem_image.2 ⟨x, Finset.mem_erase.2 ⟨hxne, mem_topLayer_iff.2 ⟨hxS, hxtop⟩⟩, rfl⟩
  have hx₀mem : x₀ ∈ topLayer S := mem_topLayer_iff.2 ⟨hx₀S, hx₀top⟩
  have hcard : (topLayer T).card < (topLayer S).card := by
    have h1 : (topLayer T).card ≤ (((topLayer S).erase x₀).image g).card :=
      Finset.card_le_card hsub
    have h2 : (((topLayer S).erase x₀).image g).card ≤ ((topLayer S).erase x₀).card :=
      Finset.card_image_le
    have h3 : ((topLayer S).erase x₀).card = (topLayer S).card - 1 :=
      Finset.card_erase_of_mem hx₀mem
    have h4 : 1 ≤ (topLayer S).card := Finset.card_pos.2 ⟨x₀, hx₀mem⟩
    omega
  rw [dSum_eq, dSum_eq, heq]
  exact Nat.mul_lt_mul_of_lt_of_le hcard (le_refl _) hm0

/-! ## ★★★★★★★★★★入れ子帰納法そのもの

原文 (NCBelyi p.5):
> on m(S), and, for each fixed value of m(S), we induct on d(S). If m(S), d(S)

★原文の『`m(S)` について帰納し、`m(S)` を固定するごとに `d(S)` について帰納する』を
**辞書式の整礎帰納法**として取る。
-/

/-- ★★★★★★★★★**入れ子帰納法の骨格**。

原文 (NCBelyi p.5):
> on m(S), and, for each fixed value of m(S), we induct on d(S). If m(S), d(S)

`m(S)` について外側、`d(S)` について内側の 2 重帰納である。
★`ℕ × ℕ` の辞書式順序を持ち出さず、**`ℕ` の帰納を 2 段重ねて**取った
——上界 `n`, `k` を introduce するだけで済む。 -/
theorem nested_induction {P : Finset ℂ → Prop}
    (base : ∀ S : Finset ℂ, maxRedDeg S = 0 → P S)
    (step : ∀ S : Finset ℂ, 0 < maxRedDeg S →
      (∀ T : Finset ℂ, maxRedDeg T < maxRedDeg S → P T) →
      (∀ T : Finset ℂ, maxRedDeg T = maxRedDeg S → dSum T < dSum S → P T) → P S) :
    ∀ S : Finset ℂ, P S := by
  have key : ∀ n k : ℕ, ∀ S : Finset ℂ, maxRedDeg S ≤ n → dSum S ≤ k → P S := by
    intro n
    induction n with
    | zero => intro k S hm _; exact base S (by omega)
    | succ n ih =>
      intro k
      induction k with
      | zero =>
        intro S hm hd
        by_cases h0 : maxRedDeg S = 0
        · exact base S h0
        · exact step S (by omega) (fun T hT => ih (dSum T) T (by omega) le_rfl)
            (fun T _ hTd => absurd hTd (by omega))
      | succ k ihk =>
        intro S hm hd
        by_cases h0 : maxRedDeg S = 0
        · exact base S h0
        · exact step S (by omega) (fun T hT => ih (dSum T) T (by omega) le_rfl)
            (fun T hTm hTd => ihk T (by omega) (by omega))
  intro S
  exact key (maxRedDeg S) (dSum S) S le_rfl le_rfl

/-- ★★★★★★★★★★**[NCBelyi] Lemma 2.4 の帰納の本体**。

原文 (NCBelyi p.5):
> hypothesis, and composing the resulting morphisms yields

`m(S) > 0` のとき **1 段下げる操作**(`g` と余分な点 `E`、落ちる点 `x₀`)が作れて、
その結論が `S` へ戻せるなら、**すべての `S` で結論が成り立つ**。

★これが原文の『replacing S by S′, β by f₀(β), applying the induction hypothesis,
and composing the resulting morphisms』の骨格である。
★★測度が下がることは `measure_lt`(第 408)が保証する。

★★★**残っているのは `descend` を実際に作ること**である
——`f₀` を `α₀` の最小多項式に取り、`E ≝ f₀(S₀)` とする段。
その定量的な部分(`f₀(β) ∉ S′`)は `CoeffBound.lean` に、
次数の降下は `ReducedDegree.lean` にある。 -/
theorem nested_induction_descend {P : Finset ℂ → Prop}
    (base : ∀ S : Finset ℂ, maxRedDeg S = 0 → P S)
    (descend : ∀ S : Finset ℂ, 0 < maxRedDeg S →
      ∃ (E : Finset ℂ) (g : ℂ → ℂ) (x₀ : ℂ),
        (∀ x ∈ S, redDeg (g x) ≤ redDeg x)
        ∧ (∀ y ∈ E, redDeg y < maxRedDeg S)
        ∧ x₀ ∈ S ∧ redDeg x₀ = maxRedDeg S ∧ redDeg (g x₀) < maxRedDeg S
        ∧ (P (S.image g ∪ E) → P S)) :
    ∀ S : Finset ℂ, P S := by
  classical
  refine nested_induction base (fun S hS ihm ihd => ?_)
  obtain ⟨E, g, x₀, hno, hE, hx₀S, hx₀top, hx₀drop, himp⟩ := descend S hS
  refine himp ?_
  rcases measure_lt S E g hno hE hx₀S hx₀top hx₀drop hS with h | ⟨h1, h2⟩
  · exact ihm _ h
  · exact ihd _ h1 h2

/-- ★★★**基底段の中身** —— `m(S) = 0` は「`S` がすべて有理点」である。

原文 (NCBelyi p.5):
> f and only if m(S) = 0 if and only if S ⊆P1(Q).

★帰納法が止まるところがちょうど `Lemma 2.2` の適用条件である。 -/
theorem forall_isRat_of_maxRedDeg_eq_zero {S : Finset ℂ} (hS : maxRedDeg S = 0)
    (hint : ∀ x ∈ S, IsIntegral ℚ x) :
    ∀ x ∈ S, ∃ q : ℚ, x = (q : ℂ) := by
  intro x hx
  exact (redDeg_eq_zero_iff (hint x hx)).1 (maxRedDeg_eq_zero_iff.1 hS x hx)

/-! ## ★出典の紐付け(`.src`) -/

def maxRedDeg.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Lemma 2.4(m(S) —— reduced degree の最大値)",
    sectionId := "ncbelyi-lemma-2-4" }

def dSum.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Lemma 2.4(d(S) —— 最大値を取る点の reduced degree の和)",
    sectionId := "ncbelyi-lemma-2-4" }

def dSum_eq_zero_iff_of_nonempty.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Lemma 2.4(d(S) = 0 ⟺ m(S) = 0 ⟺ S ⊆ ℙ¹(ℚ))",
    sectionId := "ncbelyi-lemma-2-4" }

def measure_lt.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Lemma 2.4(nested induction —— m(S′) < m(S) または m(S′) = m(S), d(S′) < d(S))",
    sectionId := "ncbelyi-lemma-2-4" }

def nested_induction.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Lemma 2.4(nested induction —— m(S) の外側・d(S) の内側の 2 重帰納)",
    sectionId := "ncbelyi-lemma-2-4" }

def nested_induction_descend.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Lemma 2.4(帰納の本体——1 段下げる操作から全体へ)",
    sectionId := "ncbelyi-lemma-2-4" }

end ABC3.Found.NCBelyi
