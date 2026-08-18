import ABC3.Found.Arakelov.PicSecMap

/-!
# Arakelov (B1) 第 142 ブロック —— ★★★★★★**`IsLocalizing` の `surj'`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★古典的な `f^n` 論法(Hartshorne II.5.1 / Stacks 01I3)

    ∀ y ∈ Γ(F, D f), ∃ a ∈ Γ(F,⊤), ∃ N, f^N · y = a|_{D f}

| 段 | 内容 |
|---|---|
| 1 | 各 `D(f·gᵢ)` は `Γ(F,D gᵢ)` の局所化(第 137)なので `f^{kᵢ} y|ᵢ = Aᵢ|` |
| 2 | `N := sup kᵢ` を取り `A2ᵢ := f^{N-kᵢ} Aᵢ` に揃える |
| 3 | 交わり `D(f·gᵢ·gⱼ)` へ落とすと `A2ᵢ` と `A2ⱼ` は一致する |
| 4 | ゆえに `∃ mᵢⱼ, f^{mᵢⱼ} A2ᵢ| = f^{mᵢⱼ} A2ⱼ|`(第 137 の `exists_of_eq`) |
| 5 | `M := sup mᵢⱼ` を取れば `f^M A2ᵢ` は**両立**し、層で貼れる |
| 6 | 貼った `a` について `f^{M+N} y = a|_{D f}`(被覆上で確認) |

★★**2 段の最大値**が要点である——第 1 段は「局所化の分母」、第 2 段は「貼り合わせのずれ」。

## ★★★逃げ道——`clear_value`

`set A2 := fun i => f ^ (N - kᵢ) • A i` としても、`rw [map_smul]` が
**定義を透かして** `A2 i` を展開してしまい `rw` が空回りする。

★**`clear_value A2`** で本体を消すと予定通りに動く
(等式 `hA2` は残るので必要なときに戻せる)。
-/

universe u v

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable (R : CommRingCat.{u}) (F : (Spec R).Modules)


/-- ★`D(f·(a·b)) ≤ D(f·a)`。 -/
theorem specD_mul_le' (f a b : (R : Type u)) :
    specD R (f * (a * b)) ≤ specD R (f * a) := by
  show PrimeSpectrum.basicOpen (f * (a * b)) ≤ PrimeSpectrum.basicOpen (f * a)
  rw [PrimeSpectrum.basicOpen_mul, PrimeSpectrum.basicOpen_mul, PrimeSpectrum.basicOpen_mul]
  exact le_inf inf_le_left (le_trans inf_le_right inf_le_left)

/-- ★`D(f·(a·b)) ≤ D(a) ⊓ D(b)`。 -/
theorem specD_mul_le_inf (f a b : (R : Type u)) :
    specD R (f * (a * b)) ≤ specD R a ⊓ specD R b := by
  show PrimeSpectrum.basicOpen (f * (a * b)) ≤ _
  rw [PrimeSpectrum.basicOpen_mul, PrimeSpectrum.basicOpen_mul]
  exact le_trans inf_le_right le_rfl

/-- ★`D(f·(a·b)) ≤ D(f·b)`。 -/
theorem specD_mul_le2 (f a b : (R : Type u)) :
    specD R (f * (a * b)) ≤ specD R (f * b) := by
  show PrimeSpectrum.basicOpen (f * (a * b)) ≤ PrimeSpectrum.basicOpen (f * b)
  rw [PrimeSpectrum.basicOpen_mul, PrimeSpectrum.basicOpen_mul, PrimeSpectrum.basicOpen_mul]
  exact le_inf inf_le_left (le_trans inf_le_right inf_le_right)

/-- ★`D(a·b) ≤ D(a)`。 -/
theorem specD_mul_le_left (a b : (R : Type u)) : specD R (a * b) ≤ specD R a := by
  show PrimeSpectrum.basicOpen (a * b) ≤ PrimeSpectrum.basicOpen a
  rw [PrimeSpectrum.basicOpen_mul]; exact inf_le_left

theorem exists_pow_smul_eq_res {n : ℕ} (g : Fin n → (R : Type u))
    (hspan : Ideal.span (Set.range g) = ⊤)
    (htriv : ∀ i, Nonempty ((restrictPresheafFunctor (Spec R) (specD R (g i))).obj F.val
      ≅ 𝟙_ (PresheafModulesOn (Spec R) (specD R (g i)))))
    (f : (R : Type u))
    (y : (((modulesSpecToSheaf.obj F).obj.obj (op (specD R f))) : Type u)) :
    ∃ (a : (((modulesSpecToSheaf.obj F).obj.obj (op (⊤ : (Spec R).Opens))) : Type u)) (N : ℕ),
      f ^ N • y = secMap R F ⊤ (specD R f) le_top a := by
  have step1 : ∀ i : Fin n, ∃ (ai : (((modulesSpecToSheaf.obj F).obj.obj
      (op (specD R (g i)))) : Type u)) (k : ℕ),
      f ^ k • (secMap R F (specD R f) (specD R (f * g i)) (specDle_left R f (g i)) y)
        = secMap R F (specD R (g i)) (specD R (f * g i)) (specDle R (g i) f) ai := by
    intro i
    haveI := isLocalizedModule_secRes R F (g i) f (htriv i).some
    obtain ⟨⟨ai, c⟩, hc⟩ := IsLocalizedModule.surj (S := Submonoid.powers f)
      (secMap R F (specD R (g i)) (specD R (f * g i)) (specDle R (g i) f))
      (secMap R F (specD R f) (specD R (f * g i)) (specDle_left R f (g i)) y)
    obtain ⟨k, hk⟩ := c.2
    exact ⟨ai, k, by rw [show f ^ k = (c : (R : Type u)) from hk, ← Submonoid.smul_def]; exact hc⟩
  choose A k hk using step1
  have step2 : ∀ i : Fin n,
      f ^ (Finset.univ.sup k) •
          (secMap R F (specD R f) (specD R (f * g i)) (specDle_left R f (g i)) y)
        = secMap R F (specD R (g i)) (specD R (f * g i)) (specDle R (g i) f)
          (f ^ (Finset.univ.sup k - k i) • A i) := by
    intro i
    have hle : k i ≤ Finset.univ.sup k := Finset.le_sup (Finset.mem_univ i)
    rw [map_smul, ← hk i, ← mul_smul, ← pow_add, Nat.sub_add_cancel hle]
  set N := Finset.univ.sup k with hNdef
  set A2 : ∀ i : Fin n, (((modulesSpecToSheaf.obj F).obj.obj (op (specD R (g i)))) : Type u) :=
    fun i => f ^ (N - k i) • A i with hA2
  have step3 : ∀ i j : Fin n,
      f ^ N • (secMap R F (specD R f) (specD R (f * (g i * g j)))
          (le_trans (specD_mul_le' R f (g i) (g j)) (specDle_left R f (g i))) y)
        = secMap R F (specD R (g i)) (specD R (f * (g i * g j)))
          (le_trans (specD_mul_le_inf R f (g i) (g j)) inf_le_left) (A2 i) := by
    intro i j
    have h := congrArg (secMap R F (specD R (f * g i)) (specD R (f * (g i * g j)))
      (specD_mul_le' R f (g i) (g j))) (step2 i)
    rw [map_smul, secMap_trans, secMap_trans] at h
    exact h
  have step3j : ∀ i j : Fin n,
      f ^ N • (secMap R F (specD R f) (specD R (f * (g i * g j)))
          (le_trans (specD_mul_le2 R f (g i) (g j)) (specDle_left R f (g j))) y)
        = secMap R F (specD R (g j)) (specD R (f * (g i * g j)))
          (le_trans (specD_mul_le_inf R f (g i) (g j)) inf_le_right) (A2 j) := by
    intro i j
    have h := congrArg (secMap R F (specD R (f * g j)) (specD R (f * (g i * g j)))
      (specD_mul_le2 R f (g i) (g j))) (step2 j)
    rw [map_smul, secMap_trans, secMap_trans] at h
    exact h
  have step4 : ∀ i j : Fin n,
      secMap R F (specD R (g i)) (specD R (f * (g i * g j)))
          (le_trans (specD_mul_le_inf R f (g i) (g j)) inf_le_left) (A2 i)
        = secMap R F (specD R (g j)) (specD R (f * (g i * g j)))
          (le_trans (specD_mul_le_inf R f (g i) (g j)) inf_le_right) (A2 j) :=
    fun i j => (step3 i j).symm.trans (step3j i j)
  have step5 : ∀ i j : Fin n, ∃ m : ℕ,
      f ^ m • secMap R F (specD R (g i)) (specD R (g i) ⊓ specD R (g j)) inf_le_left (A2 i)
        = f ^ m • secMap R F (specD R (g j)) (specD R (g i) ⊓ specD R (g j))
          inf_le_right (A2 j) := by
    intro i j
    haveI := isLocalizedModule_secRes_eq R F (specD R (g i) ⊓ specD R (g j))
      (specD R (f * (g i * g j))) (g i * g j) f (inf_specD R (g i) (g j)) rfl
      (specD_mul_le_inf R f (g i) (g j))
      (trivialOfLe F.val (specD_mul_le_left R (g i) (g j)) (htriv i).some)
    have heq : secMap R F (specD R (g i) ⊓ specD R (g j)) (specD R (f * (g i * g j)))
        (specD_mul_le_inf R f (g i) (g j))
        (secMap R F (specD R (g i)) (specD R (g i) ⊓ specD R (g j)) inf_le_left (A2 i))
      = secMap R F (specD R (g i) ⊓ specD R (g j)) (specD R (f * (g i * g j)))
        (specD_mul_le_inf R f (g i) (g j))
        (secMap R F (specD R (g j)) (specD R (g i) ⊓ specD R (g j)) inf_le_right (A2 j)) := by
      rw [secMap_trans, secMap_trans]
      exact step4 i j
    obtain ⟨c, hc⟩ := IsLocalizedModule.exists_of_eq (S := Submonoid.powers f)
      (f := secMap R F (specD R (g i) ⊓ specD R (g j)) (specD R (f * (g i * g j)))
        (specD_mul_le_inf R f (g i) (g j))) heq
    obtain ⟨m, hm⟩ := c.2
    exact ⟨m, by rw [show f ^ m = (c : (R : Type u)) from hm, ← Submonoid.smul_def,
      ← Submonoid.smul_def]; exact hc⟩
  choose m hm using step5
  set M := Finset.univ.sup (fun p : Fin n × Fin n => m p.1 p.2) with hMdef
  have hA2' : ∀ i : Fin n, A2 i = f ^ (N - k i) • A i := fun i => by rw [hA2]
  clear_value A2
  have hcomp : TopCat.Presheaf.IsCompatible (modulesSpecToSheaf.obj F).obj
      (fun i => specD R (g i)) (fun i => f ^ M • A2 i) := by
    intro i j
    have hle : m i j ≤ M :=
      Finset.le_sup (f := fun p : Fin n × Fin n => m p.1 p.2) (Finset.mem_univ (i, j))
    show secMap R F (specD R (g i)) _ inf_le_left (f ^ M • A2 i)
      = secMap R F (specD R (g j)) _ inf_le_right (f ^ M • A2 j)
    rw [show M = (M - m i j) + m i j from (Nat.sub_add_cancel hle).symm, pow_add, mul_smul,
      mul_smul]
    simp only [map_smul]
    rw [hm i j]
  obtain ⟨a, ha, -⟩ := TopCat.Sheaf.existsUnique_gluing' (modulesSpecToSheaf.obj F)
    (fun i => specD R (g i)) ⊤ (fun i => homOfLE le_top) (top_le_iSup_specD R g hspan)
    (fun i => f ^ M • A2 i) hcomp
  have ha' : ∀ i : Fin n, secMap R F ⊤ (specD R (g i)) le_top a = f ^ M • A2 i := ha
  refine ⟨a, M + N, ?_⟩
  refine TopCat.Sheaf.eq_of_locally_eq' (modulesSpecToSheaf.obj F)
    (fun i => specD R (f * g i)) (specD R f) (fun i => homOfLE (specDle_left R f (g i)))
    (specD_le_iSup R g hspan f) _ _ (fun i => ?_)
  show secMap R F (specD R f) (specD R (f * g i)) _ (f ^ (M + N) • y)
    = secMap R F (specD R f) (specD R (f * g i)) _ (secMap R F ⊤ (specD R f) le_top a)
  rw [map_smul, pow_add, mul_smul, step2 i, ← map_smul, ← hA2' i, ← ha' i, secMap_trans,
    secMap_trans]


/-! ## ★出典の紐付け(`.src`) -/

def exists_pow_smul_eq_res.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——IsLocalizing の surj')",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
