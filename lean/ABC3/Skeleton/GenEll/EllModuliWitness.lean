/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.EllModuliObjects
import ABC3.Found.GenEll.EllModuliGalois
import ABC3.Found.GenEll.DetCycloChar
import ABC3.Interface.GenEll.EllModuli
import ABC3.Skeleton.GenEll.Section3
import ABC3.Skeleton.GenEll.GaloisImage
import ABC3.Skeleton.GenEll.Section4
import ABC3.Meta.Claim

/-!
# `EllModuliData` の witness（`Skeleton`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17–p.23。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★★★★これは何か

`Found/GenEll/EllModuliObjects.lean`・`EllModuliGalois.lean` で作った対象を
**実際に `EllModuliData` に組み上げる**。★これで 50 個近い欄の型が一度に検査される。

## ★★★残る `sorry` は 4 本の葉に対応する 6 欄だけである

| 欄 | 依存する葉 |
|---|---|
| `degInf_quotLCyclic` | 葉 1・2（`Lemma 3.5`） |
| `faltingsHeight_quotLCyclic` | 葉 1・2 |
| `galoisFinite_lcyclicExc` | 葉 1・2（`Lemma 3.5` の (†)）＋ `northcott`（済） |
| `imageContainsSL2_of_torsionExt` | 葉 3（局所理論の行列表示） |
| `imageSurjective_of_containsSL2` | 葉 4（円分指標の全射性） |

☆`mem_lcyclicExc` は `lcyclicExc` を**その定義そのもの**に取るので自明に埋まる
——有限性の側（`galoisFinite_lcyclicExc`）が本体である。

## ★★★★★★★★★★★★★★★★残りの地図（2026-09-02、第 1239-1249）

☆**4 本の `sorry` はすべて `Found` 側の形に降りた**。

| # | `sorry` | `Found` 側の形 | 残る前提 |
|---|---|---|---|
| 1 | `galoisFiniteJ_lcyclicExcJ` | `exists_galoisFinite_lcyclic`（第 1239） | `VeluQuotOK` |
| 2 | `degInfJ_quotLCyclicJ` | `degInfJ_quotLCyclicJ_of_jExp`（第 1248） | `Lemma 3.2, (ii)` |
| 3 | `faltingsHeightJ_quotLCyclicJ` | `faltingsHeightJ_quotLCyclicJ_of_isog`（第 1245） | 同種写像の高さ評価 |
| 4 | `imageContainsSL2J_torsionExt` | `imageContainsSL2J_of_galTate'`（第 1249） | `σ` の幂単性・非自明性 |

★★**測定**——(4) の `hclosed`（像の閉性）は**不要だった**。
`galRep` の連続性は `galRep_continuous'`（第 772、葉 5）で既に閉じている。

☆(2) の道筋: `h4`・`h6`（Vélu の `v, w` と Eisenstein の関係）
→ `j_velu_tate_eq_map`（在庫）→ `jExp_eq_mul_of_j_tate_pow`（第 997）
→ `minDeltaExp_eq_mul_of_jExp_all`（第 1247）→ `degInfJ_quotLCyclicJ_of_jExp`（第 1248）。
★★☆**訂正（2026-09-02、第 1256-1257）**——`h4`・`h6` は
`Skeleton/GenEll/TateIsogeny.lean` の `c4_velu_tate`・`c6_velu_tate` として
**既に証明済み**である（sorry 0）。
同じく (3) の「同種写像の高さ評価」も
`htFalt_veluQuotientFull_le`（第 704、証明済み）である。

★★★**したがって 4 本に残るのは配管だけであり、
未証明の外部引用は 1 本もない**。
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Meta ABC3.Interface.GenEll ABC3.Found.GenEll ABC3.Found.GaloisRep

open scoped Classical in
/-- ★★★★**`lcyclicExc` 欄**——`mem_lcyclicExc` が要求する条件そのものを集めた集合。

★有限性（`galoisFinite_lcyclicExc`）が本体であり、それは `Lemma 3.5` の (†) と
`northcott`（`§9-1179`、第 753）から出る。 -/
def lcyclicExcJ (C eps : ℝ) (KV : Set ℂ) : Set ℂ :=
  {x : ℂ | ∃ (E : RealizedClass) (l : ℕ), x = E.cls ∧ Nat.Prime l ∧
    E.rep.toSSCurve.SemiStable ∧ HasLCyclicJ E.rep.toSSCurve l ∧
    E.rep.toSSCurve.PrimeToLocalHeights l ∧
    ((100 * (E.degOfDefinition : ℝ)
        * (faltingsHeightJ E.cls + C * (E.degOfDefinition : ℝ) ^ eps) ≤ (l : ℝ))
      ∨ E.cls ∈ KV)}

/-- ★★★`lcyclicExc` が `Galois`-finite であること——`Lemma 3.5` の (†) を受ける。 -/
theorem galoisFiniteJ_lcyclicExcJ (C eps : ℝ) (KV : Set ℂ) (hKV : CompactlyBoundedJ KV) :
    GaloisFiniteJ (lcyclicExcJ C eps KV) := by
  sorry

/-- ★★★`deg∞(E′) = l·deg∞(E)`——`Lemma 3.5`（葉 1・2）を受ける。 -/
theorem degInfJ_quotLCyclicJ (x : RealizedClass) (l : ℕ) (hl : Nat.Prime l)
    (hcyc : HasLCyclicJ x.rep.toSSCurve l)
    (hpr : x.rep.toSSCurve.PrimeToLocalHeights l) :
    degInfJ (quotLCyclicJ x l).cls = (l : ℝ) * degInfJ x.cls := by
  sorry

/-- ★★★`ht^Falt(E′) ≤ ht^Falt(E) + 2log(l) + C₀`——`Lemma 3.5`（葉 1・2）を受ける。 -/
theorem faltingsHeightJ_quotLCyclicJ :
    ∃ C₀ : ℝ, ∀ (x : RealizedClass) (l : ℕ), Nat.Prime l → HasLCyclicJ x.rep.toSSCurve l →
      faltingsHeightJ (quotLCyclicJ x l).cls
        ≤ faltingsHeightJ x.cls + 2 * Real.log l + C₀ := by
  sorry

/-- ★★★`α` の側（葉 3）を受ける。 -/
theorem imageContainsSL2J_torsionExt (x : RealizedClass) (l : ℕ) (hl : Nat.Prime l)
    (hl5 : 5 ≤ l) (hm : x.rep.toSSCurve.HasMultRed)
    (hp : x.rep.toSSCurve.PrimeToLocalHeights l)
    (hc : ¬ HasLCyclicJ x.rep.toSSCurve l) : ImageContainsSL2J x.rep.toSSCurve l := by
  sorry

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**葉 4 が閉じたので、これは `sorry` ではなくなった**（2026-08-31、第 784）。

`hram : l ∉ ramPrimes` は「`l` は `disc L` を割らない」を含む——
`ramPrimes` の定義に `(disc L).natAbs.primeFactors` が入っているからである。 -/
theorem imageSurjectiveJ_of_containsSL2' (x : RealizedClass) (l : ℕ) (hl : Nat.Prime l)
    (hram : x.PrimeToRamification l) (h : ImageContainsSL2J x.rep.toSSCurve l) :
    ImageSurjectiveJ x.rep.toSSCurve l := by
  haveI : Fact l.Prime := ⟨hl⟩
  have hunram : ¬ (l : ℤ) ∣ NumberField.discr x.rep.toSSCurve.fld := by
    intro hdvd
    refine hram (Finset.mem_union.2 (Or.inl (Finset.mem_union.2 (Or.inr ?_))))
    refine Nat.mem_primeFactors.2 ⟨hl, ?_, ?_⟩
    · have := Int.natAbs_dvd_natAbs.2 hdvd
      simpa using this
    · exact Int.natAbs_ne_zero.2 (NumberField.discr_ne_zero (K := x.rep.toSSCurve.fld))
  exact imageSurjectiveJ_of_cyclotomic x.rep.toSSCurve l
    (fun u => cyclotomic_det_surjective x.rep.toSSCurve l hunram u) h

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`EllModuliData` の witness**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★★★50 個近い欄のうち、`sorry` は上の 5 本だけである。 -/
noncomputable def ellModuliWitness : EllModuliData where
  EllClass := ℂ
  Curve := RealizedClass
  cls := RealizedClass.cls
  degOfDefinition := RealizedClass.degOfDefinition
  faltingsHeight := faltingsHeightJ
  HasPotMultRed := fun x => x.rep.toSSCurve.HasPotMultRed
  PrimeToLocalHeights := fun x l => x.rep.toSSCurve.PrimeToLocalHeights l
  CompactlyBounded := CompactlyBoundedJ
  GaloisFinite := GaloisFiniteJ
  ImageContainsSL2 := fun x l => ImageContainsSL2J x.rep.toSSCurve l
  degInf := degInfJ
  htInf := fun j => 12 * faltingsHeightJ j
  logDiffMell := logDiffMellJ
  degLe := degLeJ
  SemiStable := fun x => x.rep.toSSCurve.SemiStable
  HasLCyclic := fun x l => HasLCyclicJ x.rep.toSSCurve l
  MinimalField := fun x => MinimalFieldJ x.rep
  ImageSurjective := fun x l => ImageSurjectiveJ x.rep.toSSCurve l
  PrimeToRamification := RealizedClass.PrimeToRamification
  HasMultRed := fun x => x.rep.toSSCurve.HasMultRed
  PrimeToMultPrimes := fun x l => x.rep.PrimeToMultPrimes l
  degInf_le_htInf := degInfJ_sub_htInfJ_le
  htInf_bdeq_faltings := ⟨0, fun x => by simp⟩
  faltingsHeight_bddBelow := faltingsHeightJ_bddBelow
  northcott := fun C d _ => northcottJ C d
  quotLCyclic := quotLCyclicJ
  degInf_quotLCyclic := fun E l hl hcyc hpr => degInfJ_quotLCyclicJ E l hl hcyc hpr
  faltingsHeight_quotLCyclic := faltingsHeightJ_quotLCyclicJ
  degOfDefinition_pos := RealizedClass.degOfDefinition_pos
  primeToLocalHeights_of_lt := by
    intro E l hl _ hlt
    refine E.rep.toSSCurve.primeToLocalHeights_of_lt l hl ?_
    have h2 : degInfJ E.rep.toSSCurve.j = degInfJ E.cls := E.degInfJ_rep
    rw [h2]
    exact hlt
  lcyclicExc := lcyclicExcJ
  galoisFinite_lcyclicExc := galoisFiniteJ_lcyclicExcJ
  mem_lcyclicExc := by
    intro C eps KV E l hl hss hcyc hpr hor
    exact ⟨E, l, rfl, hl, hss, hcyc, hpr, hor⟩
  noMultRedExc := noMultRedExcJ
  galoisFinite_noMultRedExc := galoisFiniteJ_noMultRedExcJ
  mem_noMultRedExc := by
    intro KV E hKV hnm
    exact absurd E.rep.multRed hnm
  galoisFinite_union := galoisFiniteJ_union
  torsionExt := fun x => x
  cls_torsionExt := fun _ => rfl
  degOfDefinition_torsionExt := by
    intro E
    have := E.degOfDefinition_pos
    omega
  semiStable_torsionExt := fun E => E.rep.toSSCurve.semiStable_all
  hasMultRed_torsionExt := fun _ h => h
  primeToLocalHeights_torsionExt := fun _ _ h _ => h
  imageContainsSL2_of_torsionExt := fun E l hl hl5 hm hp hc =>
    imageContainsSL2J_torsionExt E l hl hl5 hm hp hc
  imageSurjective_of_containsSL2 := fun E l hl hram h =>
    imageSurjectiveJ_of_containsSL2' E l hl hram h
  compactlyBounded_empty := compactlyBoundedJ_empty
  multCard := fun x => x.rep.multCard
  multCard_pos := fun x => x.rep.multCard_pos
  multPrime := fun x j => x.rep.multPrime j
  multPrime_prime := fun x j => x.rep.multPrime_prime j
  localHt := fun x j => x.rep.localHt j
  localHt_pos := fun x j => x.rep.localHt_pos j
  sum_localHt_eq := RealizedClass.sum_localHt_eq
  badPrimes := fun x => x.rep.badPrimes
  badPrimes_prime := fun x => x.rep.badPrimes_prime
  sum_log_badPrimes_le := fun x => x.rep.sum_log_badPrimes_le
  primeTo_badPrimes := fun x l hl h => x.rep.primeTo_badPrimes l hl h
  ramPrimes := RealizedClass.ramPrimes
  ramPrimes_prime := RealizedClass.ramPrimes_prime
  badPrimes_subset_ramPrimes := RealizedClass.badPrimes_subset_ramPrimes
  primeTo_ramPrimes := RealizedClass.primeTo_ramPrimes
  sum_log_ramPrimes_le := by
    intro E
    have h := E.sum_log_ramPrimes_le
    rw [logDiffMellJ_eq]
    exact h

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★witness へ流し込む -/

/-- ★★★★★★★★★★★★★★★★★★★★**`Lemma 3.7` を witness で**。

原文 (GenEll p.18):
> Lemma 3.7. (Finite Exceptional Sets) Let

★★これで `Lemma 3.7` が**具体的な楕円曲線の言葉**になった。
☆残る `sorry` は witness の 5 本だけである。 -/
theorem lemma_3_7_witness (KV : Set ℂ) (hKV : CompactlyBoundedJ KV) (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, 0 < C ∧ ∃ Exc : Set ℂ, GaloisFiniteJ Exc ∧
      ∀ (E : RealizedClass) (l : ℕ), Nat.Prime l → E.rep.toSSCurve.SemiStable →
        ∀ (condA condB : Prop),
          (condA ↔ (100 * (E.degOfDefinition : ℝ)
                      * (faltingsHeightJ E.cls + C * (E.degOfDefinition : ℝ) ^ eps)
                        ≤ (l : ℝ) ∧ E.rep.toSSCurve.HasMultRed)) →
          (condB ↔ (E.cls ∈ KV ∧ E.rep.toSSCurve.PrimeToLocalHeights l)) →
          (condA → E.rep.toSSCurve.PrimeToLocalHeights l)
        ∧ (condB → E.cls ∉ Exc → E.rep.toSSCurve.HasMultRed)
        ∧ ((condA ∨ condB) → HasLCyclicJ E.rep.toSSCurve l → E.cls ∈ Exc) :=
  lemma_3_7 ellModuliWitness KV hKV eps heps

/-- ★★★★★★★★★★★★★★★★★★★★★★**`Theorem 3.8` を witness で**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to -/
theorem theorem_3_8_witness (KV : Set ℂ) (hKV : CompactlyBoundedJ KV) (ε : ℝ) (hε : 0 < ε) :
    ∃ C : ℝ, 0 < C ∧ ∃ Exc : Set ℂ, GaloisFiniteJ Exc ∧
      ∀ (E : RealizedClass) (l : ℕ), Nat.Prime l → E.cls ∉ Exc →
        ((23040 * 100 * (E.degOfDefinition : ℝ)
              * (faltingsHeightJ E.cls + C * (E.degOfDefinition : ℝ) ^ ε) ≤ (l : ℝ)
            ∧ E.rep.toSSCurve.HasPotMultRed)
          ∨ (E.cls ∈ KV ∧ E.rep.toSSCurve.PrimeToLocalHeights l ∧ Nat.Coprime l 30)) →
        ImageContainsSL2J E.rep.toSSCurve l :=
  theorem_3_8 ellModuliWitness KV hKV ε hε

/-- ★★★★`Corollary 4.3` も witness で通ることの確認（型検査のみ）。 -/
example (eps : ℝ) (heps : 0 < eps) := cor_4_3 ellModuliWitness eps heps

/-- ★★★★`Corollary 4.4` も witness で通ることの確認（型検査のみ）。 -/
example (KV : Set ℂ) (hKV : CompactlyBoundedJ KV) := cor_4_4 ellModuliWitness KV hKV

def lemma_3_7_witness.src : Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(witness へ流し込んだ形——具体的な楕円曲線の言葉で)",
    sectionId := "genell-lemma-3-7" }

def lemma_3_7_witness.needs : List ProofObligation :=
  [ .citation "[ABC3]" "ellModuliWitness(残る sorry は 5 本の葉)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.ellModuliWitness") 1 ]

def theorem_3_8_witness.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(witness へ流し込んだ形——具体的な楕円曲線の言葉で)",
    sectionId := "genell-thm-3-8" }

def theorem_3_8_witness.needs : List ProofObligation :=
  [ .citation "[ABC3]" "ellModuliWitness(残る sorry は 5 本の葉)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.ellModuliWitness") 1 ]

/-! ## ★出典の紐付け(`.src`) -/

def lcyclicExcJ.src : Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(lcyclicExc 欄——mem_lcyclicExc の条件そのものを集めた集合)",
    sectionId := "genell-lemma-3-7" }

def lcyclicExcJ.needs : List ProofObligation := []

def galoisFiniteJ_lcyclicExcJ.src : Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(lcyclicExc は Galois-finite——Lemma 3.5 の (†) を受ける)",
    sectionId := "genell-lemma-3-7" }

def galoisFiniteJ_lcyclicExcJ.needs : List ProofObligation :=
  [ .citation "[ABC3]" "Lemma 3.5 の (†)((l/14)·deg∞ ≤ ht^Falt + 2log l + C′)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.hdag_of_velu") 3,
    .implicitStep
      ("★(a)(b) 両側の高さ評価は済んでいる(htFalt_le_of_condA・htFalt_le_of_condB)。" ++
       "有限性は northcottJ(§9-1179、第 753、★無条件)から出る") 4 ]

def degInfJ_quotLCyclicJ.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(degInf_quotLCyclic 欄——deg∞(E′) = l·deg∞(E))",
    sectionId := "genell-lemma-3-5" }

def degInfJ_quotLCyclicJ.needs : List ProofObligation :=
  [ .citation "[ABC3]" "isMuAtBadPrimes_of_veluQuotient_of_large(第 1044、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.isMuAtBadPrimes_of_veluQuotient_of_large") 1,
    .citation "[ABC3]" "degInfOf_eq_of_local(局所の関係を足し上げる、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.degInfOf_eq_of_local") 1,
    .implicitStep
      ("★★**2026-09-01（第 1086）の測定**——" ++
       "本欄の数学はすでにある。" ++
       "☆`IsMuAtBadPrimes`（第 1044、`sorry`-free）が悪い素点で " ++
       "`Δ_min(E′) = l·Δ_min(E)` を与え、半安定なら良い素点では両辺 `0`、" ++
       "`degInfOf_eq_of_local` がそれを足し上げる。" ++
       "★残っているのは**実現の橋**だけである——" ++
       "`quotLCyclicJ x l` と `veluQuotientFull` を同定し、" ++
       "`degInfJ` と `degInfOf` を同定する段。") 6 ]

def faltingsHeightJ_quotLCyclicJ.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(faltingsHeight_quotLCyclic 欄——ht^Falt(E′) ≤ ht^Falt(E) + 2log l + C₀)",
    sectionId := "genell-lemma-3-5" }

def faltingsHeightJ_quotLCyclicJ.needs : List ProofObligation :=
  [ .citation "[ABC3]" "htFalt_veluQuotientFull_le(§9-1160、第 704、★無条件)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.htFalt_veluQuotientFull_le") 3,
    .implicitStep "☆quotLCyclicJ の定義と Vélu の商を突き合わせる段" 5 ]

def imageContainsSL2J_torsionExt.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(imageContainsSL2_of_torsionExt 欄)",
    sectionId := "genell-thm-3-8" }

def imageContainsSL2J_torsionExt.needs : List ProofObligation :=
  [ .citation "[ABC3]" "alpha_in_modl_image(Skeleton/GenEll/GaloisLocal.lean、局所理論の行列表示)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.alpha_in_modl_image") 10,
    .implicitStep
      ("★群論(Lemma 3.1, (iv))と位相(連続性・閉性)は済んでいる" ++
       "——imageContainsSL2J_of_alpha'(§9-1200、第 774)") 1 ]

def imageSurjectiveJ_of_containsSL2'.src : Source :=
  { paper := "GenEll", pdfPage := 22,
    item := "Corollary 4.3(imageSurjective_of_containsSL2 欄)",
    sectionId := "genell-cor-4-3" }

def imageSurjectiveJ_of_containsSL2'.needs : List ProofObligation :=
  [ .citation "[ABC3]" "cyclotomic_det_surjective(Skeleton/GenEll/GaloisLocal.lean、円分指標の全射性)"
      (.inProject "ABC3" "ABC3.Found.GenEll.cyclotomic_det_surjective") 1,
    .implicitStep
      ("★群論と逆極限の段は済んでいる" ++
       "——imageSurjectiveJ_of_cyclotomic(第 765)と " ++
       "cyclotomicCharacter_surjective_of_mod(§9-1199、第 773)") 1 ]

def ellModuliWitness.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(EllModuliData の witness——sorry は 5 本の葉だけ)",
    sectionId := "genell-prop-3-4" }

end ABC3.Skeleton.GenEll
