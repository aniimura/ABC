/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.AdicEvalNatural
import ABC3.Found.GaloisRep.TateXY

/-!
# Tate 級数の**自然性**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★`galois-equivariant-tate-uniformization` の段 2c-2

`AdicEvalNatural.lean` で `evalAdic` の自然性（段 2c-1）を取った。
★本ファイルはそれを **`tateXtail` / `tateYtail`** まで持ち上げる。

| 補題 | 内容 |
|---|---|
| `ringInverse_map` | `σ (Ring.inverse x) = Ring.inverse (σ x)`（`x` が単元なら） |
| `partialSum_map` | 部分和は `σ` と可換 |
| `adicSum_map` | ★`σ (adicSum a) = adicSum (σ ∘ a)`（`σ I ⊆ I`） |
| `tateXterm_map` | `f(t) = t/(1−t)²` の自然性 |
| `tateYterm_map` | `g(t) = t²/(1−t)³` の自然性 |
| `tateXtail_map` | ★★`∑_{m≥1} f(qᵐu)` の自然性 |
| `tateYtail_map` | ★★`∑_{m≥1} g(qᵐu)` の自然性 |

## ★★★★★★またしても「一意 ⟹ 自然」

`adicSum` も `Classical.choose` だが **`adicSum_unique` がある**。
★`AdicEvalNatural.lean` の `smodEq_map`（`σ I ⊆ I` なら合同が移る）と繋ぐだけ。

★★`Ring.inverse` は一般には環準同型と可換でないが、
**単元の上では可換**である（`1 − t` は `t ∈ I` なら単元、`isUnit_one_sub`）。

## ★残っている段（明示）

★`tateXpairE` / `tateXK` / `tatePtPair` の自然性は本ファイルの続きである
（`TateOrigin.lean` を要するので別ファイルになる）。
-/

namespace ABC3.Found.GaloisRep

/-! ## ★単元の逆元は自然 -/

/-- ★★**`Ring.inverse` は単元の上で環準同型と可換**。

★一般には成り立たない（非単元では `0` に潰れるため）。 -/
theorem ringInverse_map {R R' : Type} [CommRing R] [CommRing R'] (σ : R →+* R')
    {x : R} (hx : IsUnit x) :
    σ (Ring.inverse x) = Ring.inverse (σ x) := by
  obtain ⟨u, rfl⟩ := hx
  rw [Ring.inverse_unit]
  have hcast : σ (u : R) = ((Units.map (σ : R →* R') u : R'ˣ) : R') := rfl
  rw [hcast, Ring.inverse_unit, ← map_inv]
  rfl

/-! ## ★★adic 和の自然性 -/

/-- ★**部分和は `σ` と可換**（有限和）。 -/
theorem partialSum_map {R R' : Type} [CommRing R] [CommRing R'] (σ : R →+* R')
    (a : ℕ → R) (n : ℕ) :
    σ (partialSum a n) = partialSum (fun k => σ (a k)) n := by
  unfold partialSum
  rw [map_sum]

/-- ★★★★★★**adic 和は環準同型で自然**（`σ I ⊆ I` なら）。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★機構は `evalAdic_map` と同じ——**一意性 `adicSum_unique` 1 本**である。 -/
theorem adicSum_map {R R' : Type} [CommRing R] [CommRing R'] {I : Ideal R} {I' : Ideal R'}
    [IsAdicComplete I R] [IsAdicComplete I' R']
    (σ : R →+* R') (hσI : ∀ x ∈ I, σ x ∈ I') (a : ℕ → R) (ha : ∀ n, a n ∈ I ^ n)
    (ha' : ∀ n, σ (a n) ∈ I' ^ n) :
    adicSum (fun n => σ (a n)) ha' = σ (adicSum a ha) := by
  refine adicSum_unique _ ha' _ (fun n => ?_)
  rw [← partialSum_map σ a n]
  exact smodEq_map σ hσI (adicSum_spec a ha n)

/-! ## ★★★項の自然性 -/

/-- ★★**`f(t) = t/(1−t)²` は自然**（`t ∈ I` なら `1−t` は単元）。 -/
theorem tateXterm_map {R R' : Type} [CommRing R] [CommRing R'] {I : Ideal R}
    [IsAdicComplete I R] (σ : R →+* R') {t : R} (ht : t ∈ I) :
    σ (tateXterm t) = tateXterm (σ t) := by
  unfold tateXterm
  rw [map_mul, map_pow, ringInverse_map σ (isUnit_one_sub (I := I) ht), map_sub, map_one]

/-- ★★**`g(t) = t²/(1−t)³` は自然**。 -/
theorem tateYterm_map {R R' : Type} [CommRing R] [CommRing R'] {I : Ideal R}
    [IsAdicComplete I R] (σ : R →+* R') {t : R} (ht : t ∈ I) :
    σ (tateYterm t) = tateYterm (σ t) := by
  unfold tateYterm
  rw [map_mul, map_pow, map_pow, ringInverse_map σ (isUnit_one_sub (I := I) ht),
    map_sub, map_one]

/-! ## ★★★★★★★★尾の自然性 -/

/-- ★★★★★★★★**`∑_{m≥1} f(qᵐu)` は自然**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★`adicSum_map`（一意性）と `tateXterm_map`（項ごと）を繋ぐ。 -/
theorem tateXtail_map {R R' : Type} [CommRing R] [CommRing R'] {I : Ideal R} {I' : Ideal R'}
    [IsAdicComplete I R] [IsAdicComplete I' R']
    (σ : R →+* R') (hσI : ∀ x ∈ I, σ x ∈ I') (u q : R) (hq : q ∈ I) (hq' : σ q ∈ I') :
    σ (tateXtail u q hq) = tateXtail (σ u) (σ q) hq' := by
  have hmem : ∀ n : ℕ, q ^ (n + 1) * u ∈ I := by
    intro n
    exact Ideal.mul_mem_right u I (Ideal.pow_mem_of_mem I hq _ (Nat.succ_pos n))
  have hterm : ∀ n : ℕ, σ (tateXterm (q ^ (n + 1) * u)) = tateXterm ((σ q) ^ (n + 1) * σ u) := by
    intro n
    rw [tateXterm_map (I := I) σ (hmem n), map_mul, map_pow]
  have ha' : ∀ n : ℕ, σ (tateXterm (q ^ (n + 1) * u)) ∈ I' ^ n := by
    intro n
    rw [hterm n]
    exact tateXtail_aux hq' n
  show σ (adicSum (fun n => tateXterm (q ^ (n + 1) * u)) (tateXtail_aux hq)) = _
  rw [← adicSum_map σ hσI _ (tateXtail_aux hq) ha']
  exact adicSum_congr ha' (tateXtail_aux hq') hterm

/-- ★★★★★★★★**`∑_{m≥1} g(qᵐu)` は自然**。 -/
theorem tateYtail_map {R R' : Type} [CommRing R] [CommRing R'] {I : Ideal R} {I' : Ideal R'}
    [IsAdicComplete I R] [IsAdicComplete I' R']
    (σ : R →+* R') (hσI : ∀ x ∈ I, σ x ∈ I') (u q : R) (hq : q ∈ I) (hq' : σ q ∈ I') :
    σ (tateYtail u q hq) = tateYtail (σ u) (σ q) hq' := by
  have hmem : ∀ n : ℕ, q ^ (n + 1) * u ∈ I := by
    intro n
    exact Ideal.mul_mem_right u I (Ideal.pow_mem_of_mem I hq _ (Nat.succ_pos n))
  have hterm : ∀ n : ℕ, σ (tateYterm (q ^ (n + 1) * u)) = tateYterm ((σ q) ^ (n + 1) * σ u) := by
    intro n
    rw [tateYterm_map (I := I) σ (hmem n), map_mul, map_pow]
  have ha' : ∀ n : ℕ, σ (tateYterm (q ^ (n + 1) * u)) ∈ I' ^ n := by
    intro n
    rw [hterm n]
    exact tateYtail_aux hq' n
  show σ (adicSum (fun n => tateYterm (q ^ (n + 1) * u)) (tateYtail_aux hq)) = _
  rw [← adicSum_map σ hσI _ (tateYtail_aux hq) ha']
  exact adicSum_congr ha' (tateYtail_aux hq') hterm

/-! ### ★出典の紐付け(`.src`)

★★`Definition 3.3`（Tate 一意化）の**段 2c-2 の前半**である。 -/

def adicSum_map.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——adic 和の自然性)",
    sectionId := "genell-def-3-3" }

def tateXtail_map.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——X の尾の自然性)",
    sectionId := "genell-def-3-3" }

def tateYtail_map.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——Y の尾の自然性)",
    sectionId := "genell-def-3-3" }

def tateXtail_map.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "adicSum_unique(adic 和の一意性)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.adicSum_unique") 15,
    .citation "[ABC3]" "smodEq_map(σ I ⊆ I なら合同が移る)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.smodEq_map") 15,
    .citation "[ABC3]" "isUnit_one_sub(t ∈ I なら 1−t は単元)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.isUnit_one_sub") 15,
    .implicitStep
      ("★★★Ring.inverse は一般には環準同型と可換でないが、" ++
       "**単元の上では可換**である。★1−t は t ∈ I なら単元なので使える") 15,
    .implicitStep
      ("★★残る段: tateXpairE / tateXK / tatePtPair の自然性" ++
       "(TateOrigin.lean を要するので別ファイル)") 15 ]

end ABC3.Found.GaloisRep
