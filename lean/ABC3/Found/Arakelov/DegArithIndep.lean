/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.DegArith
import ABC3.Meta.Claim

/-!
# 段 D が閉じた —— **`deg_F` は切断の取り方に依らない**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.4。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

## ★★★★★★★★★★これは何か

`§9-776` が与えたのは「**大域函数倍で不変**」であった。
本ファイルは残りの段を埋め、**任意の二つの非零切断が同じ次数を与える**ことを示す。

## ★★★★★機構（3 段）

    (1) norm_congr_of_sheafify_eq : アルキメデスのノルムは**層化した切断だけで決まる**
    (2) exists_smul_mem_span_gamma: 層化した側では階数 1 なので a·s と r·t が一致する（段 B）
    (3) degArith_smul             : 大域函数倍での不変性（§9-776）

★(1) が鍵である。前層の水準では `Γ(L)` に「局所的に 0」な切断が残りうるが、
**そのような切断はどの複素点でもノルム 0 を与える**——
`AMetric.norm` は自明化を通した値を点で評価するだけで、
`norm_eq`（チャート独立性）によって**十分小さい開集合で測ってよい**からである。
★★したがって `L.norm s p` は `s` の層化した像だけで決まる。
★★★これで「前層の水準の逸脱」が**この段では無害である**ことが示された。

★(2) は段 B の `exists_smul_mem_span_gamma`（層化した `Γ` の階数 1 性）である。
非零因子 `r` と `a` を得て、`a·η(s) = r·η(t)` とする。
`a ≠ 0` は在庫の `smul_ne_zero_gamma`（層化した `Γ` は捻れを持たない）から出る
——`Γ` は可逆（`invertible_gamma_toInvSheaf`）ゆえ平坦だからである。

## ★★★★★★★★★★結論

    `degArith F L s = degArith F L t`   （`s`, `t` は層化した像が非零、ノルムはどの埋め込みでも非零）

★★★これで台帳 `arakelov-degF-finite-places` の**段 D が閉じた**。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite NumberField
open ABC3.Found.GenEll

/-! ## ★★★★★★(1) ノルムは層化した切断だけで決まる -/

/-- ★`⊤` からの制限は経由する開集合を通しても同じ（前層の関手性）。 -/
theorem secOn_le {X : Scheme.{0}} (M : X.PresheafOfModules) {V W : X.Opens} (hVW : V ≤ W)
    (s : (M.obj (op ⊤) : Type)) :
    secOn M V s = (M.map (homOfLE hVW).op).hom (secOn M W s) := by
  show (M.map (homOfLE (le_top)).op).hom s = _
  exact congrArg (fun (m : _ ⟶ _) => (ModuleCat.Hom.hom m) s)
    (M.map_comp (homOfLE (le_top : W ≤ ⊤)).op (homOfLE hVW).op)

/-- ★★★★**ノルムは切断の「その開集合への制限」だけで決まる**。

★`trivValue M V e s = trivEquiv M V e (secOn M V s)` であり、
`secOn` しか見ていないからである。 -/
theorem norm_congr_of_secOn_eq {X : Scheme.{0}} (L : AMetric X)
    (s t : (L.sheaf.obj (op ⊤) : Type)) (p : Spec (CommRingCat.of ℂ) ⟶ X)
    (V : X.Opens) (hp : p ⁻¹ᵁ V = ⊤)
    (e : (restrictPresheafFunctor X V).obj L.sheaf ≅ 𝟙_ (PresheafModulesOn X V))
    (h : secOn L.sheaf V s = secOn L.sheaf V t) :
    L.norm s p = L.norm t p := by
  rw [AMetric.norm_eq L s p V e hp, AMetric.norm_eq L t p V e hp]
  have hv : trivValue L.sheaf V e s = trivValue L.sheaf V e t := by
    rw [trivValue_eq_trivEquiv, trivValue_eq_trivEquiv, h]
  show ‖evalOn p V hp (trivValue L.sheaf V e s)‖ * _
    = ‖evalOn p V hp (trivValue L.sheaf V e t)‖ * _
  rw [hv]

/-- ★★★★★★★★★★**アルキメデスのノルムは層化した切断だけで決まる**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★★★これが「前層の水準で作業していること」がアルキメデス側で**無害である**ことの証明である。
機構は次の 2 つ:
  ・層化の単位は**局所単射**である（`isLocallyInjective_sheafifyUnit`）——
    像が等しければ、ある被覆の上で切断そのものが等しい。
  ・`AMetric.norm` は**チャート独立**である（`norm_eq`）——
    その被覆の元とチャートの共通部分で測ってよい。 -/
theorem norm_congr_of_sheafify_eq {X : Scheme.{0}} (L : AMetric X)
    (s t : (L.sheaf.obj (op ⊤) : Type))
    (h : ((sheafifyUnit X L.sheaf).app (op ⊤)).hom s
       = ((sheafifyUnit X L.sheaf).app (op ⊤)).hom t)
    (p : Spec (CommRingCat.of ℂ) ⟶ X) :
    L.norm s p = L.norm t p := by
  have hmem := @Presheaf.equalizerSieve_mem _ _ _ _ _ _ _ _ (Opens.grothendieckTopology X) _ _
    ((PresheafOfModules.toPresheaf X.ringCatSheaf.obj).map (sheafifyUnit X L.sheaf))
    (isLocallyInjective_sheafifyUnit X L.sheaf) (op (⊤ : X.Opens)) s t h
  obtain ⟨W, i, hi, hWmem⟩ := hmem (p.base default) trivial
  obtain ⟨c⟩ := nonempty_normChart L.triv p
  have hWtop : p ⁻¹ᵁ W = ⊤ :=
    preimage_eq_top_of_mem p W (fun z => by rw [Subsingleton.elim z default]; exact hWmem)
  have hpV : p ⁻¹ᵁ (W ⊓ c.V) = ⊤ := by
    show p ⁻¹ᵁ W ⊓ p ⁻¹ᵁ c.V = ⊤
    rw [hWtop, c.hp, inf_idem]
  have hsecW : secOn L.sheaf W s = secOn L.sheaf W t := by
    have hi' : (L.sheaf.map (homOfLE (le_top : W ≤ ⊤)).op).hom s
        = (L.sheaf.map (homOfLE (le_top : W ≤ ⊤)).op).hom t := by
      have h2 := hi
      rwa [Subsingleton.elim i (homOfLE (le_top : W ≤ ⊤))] at h2
    exact hi'
  have hsec : secOn L.sheaf (W ⊓ c.V) s = secOn L.sheaf (W ⊓ c.V) t := by
    rw [secOn_le L.sheaf (inf_le_left : W ⊓ c.V ≤ W) s,
      secOn_le L.sheaf (inf_le_left : W ⊓ c.V ≤ W) t, hsecW]
  exact norm_congr_of_secOn_eq L s t p (W ⊓ c.V) hpV
    (trivialOfLe L.sheaf (inf_le_right : W ⊓ c.V ≤ c.V) c.e) hsec

/-- ★★**アルキメデス次数は点ごとのノルムだけで決まる**。 -/
theorem archDeg_congr (F : Type) [Field F] [NumberField F]
    (L : AMetric (Spec (CommRingCat.of (𝓞 F))))
    (s t : (L.sheaf.obj (op ⊤) : Type))
    (h : ∀ p : Spec (CommRingCat.of ℂ) ⟶ Spec (CommRingCat.of (𝓞 F)),
      L.norm s p = L.norm t p) :
    archDeg F L s = archDeg F L t := by
  show -(∑ σ : F →+* ℂ, Real.log (L.norm s (embSpecPoint F σ))) / _
      = -(∑ σ : F →+* ℂ, Real.log (L.norm t (embSpecPoint F σ))) / _
  simp only [h]

/-! ## ★★★★★★★(2) 層化した側の捧れの無さ（在庫） -/

/-- ★`Γ-Spec` 同型は全射である（前層の側に持ち上げられる）。 -/
theorem exists_preimage_gammaSpec (R : CommRingCat.{0}) (r : (R : Type)) :
    ∃ c : (Γ(Spec R, (⊤ : (Spec R).Opens)) : Type),
      (Scheme.ΓSpecIso R).hom.hom c = r :=
  ⟨(Scheme.ΓSpecIso R).inv.hom r,
    congrArg (fun (m : _ ⟶ _) => CommRingCat.Hom.hom m r) (Scheme.ΓSpecIso R).inv_hom_id⟩

/-! ## ★★★★★★★★★★(3) 段 D —— 切断の取り方に依らない -/

/-- ★★★★★★★★★★**段 D が閉じた**——`deg_F` は切断の取り方に依らない。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★機構:
  ・段 B の `exists_smul_mem_span_gamma` で `a·η(s) = r·η(t)`（`r` は非零因子）を得る。
  ・`a ≠ 0` は平坦性（捻れの無さ）から出る。
  ・`Γ-Spec` 同型で `a`, `r` を前層の側に持ち上げ、`ca·s` と `cr·t` を作る。
  ・**層化した像が等しい**ので、有限側は自明に等しく、
    アルキメデス側も `norm_congr_of_sheafify_eq` で等しい。
  ・最後に `degArith_smul`（§9-776）で `ca`, `cr` を落とす。

★★★これで台帳 `arakelov-degF-finite-places` の**段 D が閉じた**。 -/
theorem degArith_congr (F : Type) [Field F] [NumberField F]
    (L : AInv (Spec (CommRingCat.of (𝓞 F))))
    (s t : (L.carrier.sheaf.obj (op ⊤) : Type))
    (hgs : gammaSheafifyM (CommRingCat.of (𝓞 F)) L s ≠ 0)
    (hgt : gammaSheafifyM (CommRingCat.of (𝓞 F)) L t ≠ 0)
    (hs : ∀ σ : F →+* ℂ, L.carrier.norm s (embSpecPoint F σ) ≠ 0)
    (ht : ∀ σ : F →+* ℂ, L.carrier.norm t (embSpecPoint F σ) ≠ 0) :
    degArith F L s = degArith F L t := by
  obtain ⟨r, hr, hmem⟩ := exists_smul_mem_span_gamma (CommRingCat.of (𝓞 F)) L
      (gammaSheafifyM (CommRingCat.of (𝓞 F)) L s)
      (gammaSheafifyM (CommRingCat.of (𝓞 F)) L t) hgs
  obtain ⟨a, ha⟩ := Submodule.mem_span_singleton.mp hmem
  have hr0 : r ≠ 0 := nonZeroDivisors.ne_zero hr
  have hrt : r • gammaSheafifyM (CommRingCat.of (𝓞 F)) L t ≠ 0 :=
    smul_ne_zero_gamma (CommRingCat.of (𝓞 F)) L _ hgt r hr0
  have ha0 : a ≠ 0 := by
    intro h0
    rw [h0, zero_smul] at ha
    exact hrt ha.symm
  obtain ⟨ca, hca⟩ := exists_preimage_gammaSpec (CommRingCat.of (𝓞 F)) a
  obtain ⟨cr, hcr⟩ := exists_preimage_gammaSpec (CommRingCat.of (𝓞 F)) r
  have key : gammaSheafifyM (CommRingCat.of (𝓞 F)) L (ca • s)
      = gammaSheafifyM (CommRingCat.of (𝓞 F)) L (cr • t) := by
    rw [gammaSheafifyM_smul, gammaSheafifyM_smul, hca, hcr, ha]
  have hnorm : ∀ p : Spec (CommRingCat.of ℂ) ⟶ Spec (CommRingCat.of (𝓞 F)),
      L.carrier.norm (ca • s) p = L.carrier.norm (cr • t) p :=
    fun p => norm_congr_of_sheafify_eq L.carrier (ca • s) (cr • t) key p
  have hA : archDeg F L.carrier (ca • s) = archDeg F L.carrier (cr • t) :=
    archDeg_congr F L.carrier (ca • s) (cr • t) hnorm
  have hmain : degArith F L (ca • s) = degArith F L (cr • t) := by
    show degFin (CommRingCat.of (𝓞 F)) L
        (gammaSheafifyM (CommRingCat.of (𝓞 F)) L (ca • s)) / _
      + archDeg F L.carrier (ca • s) = _
    rw [key, hA]
    rfl
  rw [← degArith_smul F L ca s (by rw [hca]; exact ha0) hgs hs,
      ← degArith_smul F L cr t (by rw [hcr]; exact hr0) hgt ht, hmain]

/-! ### ★出典の紐付け(`.src`) -/

def norm_congr_of_sheafify_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(アルキメデスのノルムは層化した切断だけで決まる)",
    sectionId := "genell-def-1-1-ii" }

def degArith_congr.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(段 D が閉じた——deg_F は切断の取り方に依らない)",
    sectionId := "genell-def-1-1-ii" }

def degArith_congr.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "degArith_smul(大域函数倍での不変性、§9-776)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.degArith_smul") 4,
    .citation "[ABC3]" "exists_smul_mem_span_gamma(層化した Γ の階数 1 性、段 B)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.exists_smul_mem_span_gamma") 4,
    .citation "[ABC3]" "isLocallyInjective_sheafifyUnit(層化の単位は局所単射)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.isLocallyInjective_sheafifyUnit") 4,
    .citation "[ABC3]" "AMetric.norm_eq(ノルムはチャートの取り方に依らない)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.AMetric.norm_eq") 4,
    .citation "[ABC3]" "smul_ne_zero_gamma(層化した Γ は捻れを持たない、在庫)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.smul_ne_zero_gamma") 4,
    .implicitStep
      ("★原文は deg_F が切断の取り方に依らないことを明示しない" ++
       "——[Szp] Prop 1.1 の同型を引くことで済ませている。" ++
       "★★本ブロックはその同型を経由せず、古典的な定義から直接示した") 4 ]

end ABC3.Found.Arakelov
