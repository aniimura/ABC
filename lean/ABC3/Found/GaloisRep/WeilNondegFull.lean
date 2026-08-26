import ABC3.Found.GaloisRep.DegBoundTower
import ABC3.Found.GaloisRep.NondegStep

/-!
# Galois (G5) 第 327 ブロック —— **★★★★★★★★★Weil 対の非退化性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★到達点

> **`e_n(S, ·) ≡ 1 ⟹ S = 0`**(`weilPairing_nondegenerate`)

★★★これが (G5) に残っていた**最後の数学**である。第 197 で「非退化性は
`F(E)^{E[n]} = [n]^*F(E)` の 1 つに絞れた」と測って以来の穴が埋まる。

## ★★★★★★★`hfix` の証明——5 つの部品の配線

    L := F(x_n, y_n)
    (a) x は L 上整、最小多項式の次数 ≤ n²   第 325 `finrank_adjoin_le_of_monic_rel`
    (b) y ∈ L(x)                             第 325 `coordY_mem_of_mem`
    (c) ゆえに L(x) = F(E)                   第 326 `intermediateField_eq_top_of_coord`
    (d) L ⊆ Fix、[F(E):Fix] = n² で挟む      第 326 `eq_of_finrank_bound` + 第 196 Artin
    (e) L ⊆ μF(E)                            第 118 `pointHom_genX/genY`

★`Fix` は `FixedPoints.subfield (torsGroup W n) F(E)`(第 196)を
`Subfield.toIntermediateField` で中間体に持ち上げた
——`F` の元が固定されるのは、群が **`F` 代数自己同型**で作用するからである。

## ★★★★★★当初の道との違い

第 196 は「`deg[n] = n²` には**分点多項式**か**双対同種**のどちらかが要る
——本プロジェクトがこれまで避けてきた量である」と書いた。
★★実際に効いたのは**分点多項式の側**だが、`x([n]P) = Φ_n/ΨSq_n` の
**帰納は第 42-52 で既に回してあった**(`MulOK`)。
★★★★★避けていたのは量ではなく**その帰納**であり、それは在庫にあった。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `hfix_of_mulByN` | ★★★★★★★★**`F(E)^{E[n]} = [n]^*F(E)`** |
| `exists_pairing_ne_one` | ★★★★★★★`S ≠ 0` なら対が `1` でない `Q` がある |
| `weilPairing_nondegenerate` | ★★★★★★★★★**非退化性** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IntermediateField

variable {F : Type} [Field F] [DecidableEq F] [IsAlgClosed F] [Infinite F] [CharZero F]
  (W : WeierstrassCurve.Affine F) [W.IsElliptic] [IsDedekindDomain W.CoordinateRing]

/-! ## ★★★★★★★★`F(E)^{E[n]} = [n]^*F(E)` -/

set_option maxHeartbeats 3200000 in
/-- ★★★★★★★★**`E[n]` で固定される関数は `[n]` の引き戻しである**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`L := F(x_n, y_n)` を取り、`[F(E) : L] ≤ n²` と Artin の `[F(E) : Fix] = n²` で挟む。 -/
theorem hfix_of_mulByN (h4 : (4 : F) ≠ 0) (hΔ : W.Δ ≠ 0) (n : ℕ) (hn : 1 ≤ n)
    {μ : W.CoordinateRing →+* W.FunctionField} (hinj : Function.Injective μ)
    (hμF : ∀ d : F, μ (algebraMap F W.CoordinateRing d) = algebraMap F W.FunctionField d)
    {xn yn : W.FunctionField} (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular xn yn)
    (hμP : n • genericPoint W = Point.some xn yn hns)
    (hμx : μ (genX W) = xn) (hμy : μ (genY W) = yn)
    (hxeq : xn * Polynomial.eval₂ (algebraMap F W.FunctionField) (coordX W) (W.ΨSq (n : ℤ))
        = Polynomial.eval₂ (algebraMap F W.FunctionField) (coordX W) (W.Φ (n : ℤ)))
    (hyeq : (yn - (W.map (algebraMap F W.FunctionField)).negY xn yn)
          * Polynomial.eval₂ (algebraMap F W.FunctionField) (coordX W) (W.ΨSq (n : ℤ)) ^ 2
        = Polynomial.eval₂ (algebraMap F W.FunctionField) (coordX W) (W.preΨ (2 * (n : ℤ)))
          * (coordY W - (W.map (algebraMap F W.FunctionField)).negY (coordX W) (coordY W))) :
    ∀ z : W.FunctionField,
      (∀ Q : W.Point, n • Q = 0 → translateAut W h4 Q z = z) → ∃ v, z = muExt W hinj v := by
  letI := MulSemiringAction.compHom W.FunctionField (torsionAutHom W h4 n)
  set L : IntermediateField F W.FunctionField := IntermediateField.adjoin F {xn, yn} with hLdef
  have hxnL : xn ∈ L := IntermediateField.subset_adjoin F _ (by simp)
  have hynL : yn ∈ L := IntermediateField.subset_adjoin F _ (by simp)
  have hdlt : (W.ΨSq (n : ℤ)).natDegree < (W.Φ (n : ℤ)).natDegree := by
    have h1 : (W.Φ (n : ℤ)).natDegree = n ^ 2 := by rw [natDegree_Φ_eq]; simp
    have h2 : (W.ΨSq (n : ℤ)).natDegree ≤ n ^ 2 - 1 := by
      have h := W.natDegree_ΨSq_le (n : ℤ); simpa using h
    have h3 : 1 ≤ n ^ 2 := Nat.one_le_pow 2 n (by omega)
    omega
  obtain ⟨hint, hdeg⟩ := finrank_adjoin_le_of_monic_rel L (coordX W) xn hxnL
    (W.Φ (n : ℤ)) (W.ΨSq (n : ℤ)) (monic_Φ W _) hdlt hxeq
  have hdeg2 : (minpoly L (coordX W)).natDegree ≤ n ^ 2 := by
    have h1 : (W.Φ (n : ℤ)).natDegree = n ^ 2 := by rw [natDegree_Φ_eq]; simp
    omega
  have hxmem : coordX W ∈ IntermediateField.restrictScalars F (L⟮coordX W⟯) :=
    IntermediateField.mem_adjoin_simple_self _ _
  have hLsub : ∀ u : W.FunctionField, u ∈ L →
      u ∈ IntermediateField.restrictScalars F (L⟮coordX W⟯) := by
    intro u hu
    exact (L⟮coordX W⟯).algebraMap_mem ⟨u, hu⟩
  have hymem : coordY W ∈ IntermediateField.restrictScalars F (L⟮coordX W⟯) :=
    coordY_mem_of_mem W _ n hn hyeq (hLsub xn hxnL) (hLsub yn hynL) hxmem
  have htop : IntermediateField.restrictScalars F (L⟮coordX W⟯) = ⊤ :=
    intermediateField_eq_top_of_coord W _ hxmem hymem
  obtain ⟨hle, hpos⟩ := finrank_le_of_adjoin_top L (coordX W) hint (n ^ 2) hdeg2 htop
  have hFmem : ∀ c : F, algebraMap F W.FunctionField c ∈
      FixedPoints.subfield (torsGroup W n) W.FunctionField := by
    intro c g
    show (torsionAutHom W h4 n g) (algebraMap F W.FunctionField c) = _
    exact AlgEquiv.commutes _ c
  set Fx : IntermediateField F W.FunctionField :=
    (FixedPoints.subfield (torsGroup W n) W.FunctionField).toIntermediateField hFmem with hFxdef
  have hmuFix : ∀ w : W.FunctionField, muExt W hinj w ∈ Fx := fun w =>
    muExt_mem_fixedPoints W h4 n hinj hμF hns hμx hμy hμP w
  have hLFx : L ≤ Fx := by
    rw [hLdef]
    refine IntermediateField.adjoin_le_iff.2 ?_
    rintro u hu
    simp only [Set.mem_insert_iff, Set.mem_singleton_iff] at hu
    rcases hu with rfl | rfl
    · have h := hmuFix (coordX W); rwa [muExt_coordX, hμx] at h
    · have h := hmuFix (coordY W); rwa [muExt_coordY, hμy] at h
  have hFxrank : Module.finrank Fx W.FunctionField = n ^ 2 :=
    finrank_fixedPoints_torsGroup W h4 hΔ n hn (fun k hk1 _ => Nat.cast_ne_zero.2 (by omega))
  have hLeq : L = Fx := eq_of_finrank_bound L Fx hLFx (n ^ 2) hpos hle hFxrank
  have hmuF2 : ∀ c : F, algebraMap F W.FunctionField c ∈ (muExt W hinj).fieldRange := by
    intro c
    refine ⟨algebraMap F W.FunctionField c, ?_⟩
    have h := muExt_algebraMap W hinj (algebraMap F W.CoordinateRing c)
    rw [hμF c] at h
    rwa [← IsScalarTower.algebraMap_apply] at h
  set Mu : IntermediateField F W.FunctionField :=
    (muExt W hinj).fieldRange.toIntermediateField hmuF2 with hMudef
  have hLMu : L ≤ Mu := by
    rw [hLdef]
    refine IntermediateField.adjoin_le_iff.2 ?_
    rintro u hu
    simp only [Set.mem_insert_iff, Set.mem_singleton_iff] at hu
    rcases hu with rfl | rfl
    · show u ∈ Mu
      rw [hMudef]
      exact ⟨coordX W, by rw [muExt_coordX, hμx]⟩
    · show u ∈ Mu
      rw [hMudef]
      exact ⟨coordY W, by rw [muExt_coordY, hμy]⟩
  intro z hz
  have hzFx : z ∈ Fx := by
    rw [hFxdef]
    intro g
    show translateAut W h4 _ z = z
    exact hz _ ((mem_ker_nsmulEndo n _).1 (Multiplicative.toAdd g).2)
  rw [← hLeq] at hzFx
  have hzMu := hLMu hzFx
  rw [hMudef] at hzMu
  obtain ⟨v, hv⟩ := hzMu
  exact ⟨v, hv.symm⟩

/-! ## ★★★★★★★★★非退化性 -/

set_option maxHeartbeats 3200000 in
/-- ★★★★★★★`S ≠ 0` なら `e_n(S, Q) ≠ 1` となる `n` 等分点 `Q` がある。 -/
theorem exists_pairing_ne_one (h2 : IsUnit (2 : F)) (h4 : (4 : F) ≠ 0) (n : ℕ) (hn : 1 ≤ n)
    (x y : F) (hPx : W.Nonsingular x y) (hzero : n • Point.some x y hPx = 0)
    (hne : Point.some x y hPx ≠ 0) :
    ∃ Q, n • Q = 0 ∧ weilPairingVal W n (Point.some x y hPx) Q ≠ 1 := by
  classical
  obtain ⟨x', y', h', heq, hx, hy⟩ :=
    mulOK_of_ne (W.map (algebraMap F W.FunctionField)) (nonsingular_coord W) n hn
      (fun k hk1 _ => psiSq_coordX_ne_zero W k hk1)
  rw [WeierstrassCurve.map_ΨSq, Polynomial.eval_map, WeierstrassCurve.map_Φ,
    Polynomial.eval_map] at hx
  rw [WeierstrassCurve.map_ΨSq, Polynomial.eval_map, WeierstrassCurve.map_preΨ,
    Polynomial.eval_map] at hy
  have htr : Transcendental F x' := transcendental_mulByN_coordX W n hn hx
  have hinj : Function.Injective (pointHom W h'.1) :=
    pointHom_injective_of_transcendental W h'.1 (fun p hp hp0 => htr ⟨p, hp, hp0⟩)
  have hfix := hfix_of_mulByN W h4 W.isUnit_Δ.ne_zero n hn hinj
    (pointHom_algebraMap W h'.1) h' heq (pointHom_genX W h'.1) (pointHom_genY W h'.1) hx hy
  exact nondeg_of_fixedField W h2 h4 n hn (fun k hk1 _ => Nat.cast_ne_zero.2 (by omega))
    hinj (pointHom_algebraMap W h'.1) h' heq (pointHom_genX W h'.1) (pointHom_genY W h'.1)
    hfix hPx hzero hne

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★**Weil 対の非退化性**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★これが (G5) に残っていた最後の数学である。 -/
theorem weilPairing_nondegenerate (h2 : IsUnit (2 : F)) (h4 : (4 : F) ≠ 0) (n : ℕ) (hn : 1 ≤ n)
    (S : W.Point) (hS : n • S = 0)
    (hall : ∀ R : W.Point, n • R = 0 → weilPairingVal W n S R = 1) : S = 0 := by
  by_contra hne
  match S, hS, hne, hall with
  | Point.some x y hPx, hS, hne, hall =>
    obtain ⟨Q, hQ, hQne⟩ := exists_pairing_ne_one W h2 h4 n hn x y hPx hS hne
    exact hQne (hall Q hQ)

/-! ## ★出典の紐付け(`.src`) -/

def hfix_of_mulByN.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の非退化性——F(E)^{E[n]} = [n]^*F(E))",
    sectionId := "genell-thm-3-8" }

def weilPairing_nondegenerate.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の非退化性)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
