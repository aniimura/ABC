import ABC3.Found.PGC.NormalBasisFunctional
import Mathlib.Algebra.AlgebraicCard

/-!
# 両立汎関数 `HasCoherentFunctional` —— Maschke 上げと塔の組み立て

前段 `Found/PGC/NormalBasisFunctional.lean` は

  `SmoothModelCarrier K ⟺ HasCoherentFunctional K`

を示し、残る数学を「**1 つの** `λ : K̄ →+ K` で**全**有限次 Galois 水準を同時に良くする」
1 点に絞った。本ファイルはそのうち **(i) Maschke 上げ**を**完全に**埋め、
残りを **(ii) 可算共終塔の存在**という 1 つの仮説に還元する。

## ★★★前段の見立てより短い道(本ファイルで見つかったこと)

前段は (i) について「環論核は在る(`exists_unit_of_central_idempotent`)。
残るのは**群環 `K[G']` との同一視の配線**だけ」と見立てていた。
★**その配線は要らなかった。** 実際に必要だったのは、次の**初等的な**観察である:

`ε := ι ∘ ρ`(`ρ` は同変な収縮)が冪等で、`φ` が `B` の同変自己同型なら

  `Ψ := ι ∘ φ ∘ ρ + (1 − ι ∘ ρ)`

は `A` の同変自己同型で、逆は `φ` を `φ⁻¹` に替えたもの——★**これだけ**である
(`OrbitModel.exists_bijective_extension`)。群環も、単元の持ち上げも、中心冪等元も、
半単純性も、跡写像も出てこない。★`exists_unit_of_central_idempotent` は**使っていない**。

同様に「`K[G']` の単元の持ち上げ」に対応する部分は、
「良い汎関数どうしは同変加法同型で移り合う」(`exists_addEquiv_of_bijective`)という
Frobenius 相互律の直接の帰結に置き換わった。

## 内容

* §1 抽象核 —— 有限部分群上の**軌道和**(群と加法群だけ。体も分岐も出ない)。
* §2 抽象核 —— **Maschke 上げ**(同上)。
* §3 抽象核 —— ℕ 上の**従属選択**と**直極限**(同上)。
* §4 具体層 —— 2 つの有限水準の間の配線(★Γ_K を経由して `Algebra ↥L ↥L'` を避ける)。
* §5 具体層 —— 塔の組み立てと `HasCoherentFunctional` への還元。

## ★本ファイルが残した穴と、その後始末(2026-09-08 に同じ波で閉じた)

本ファイル自身は `HasCofinalGaloisTower K`
(「`K̄/K` の有限次 Galois 部分拡大の**増大列**で `K̄` を覆うものが取れる」)を仮説として残す。
★**その仮説は `Found/PGC/AlgebraicGenerators.lean` で証明済み**である
(`ABC3.Found.PGC.hasCofinalGaloisTower`、Krasner の補題 + `K` の可分性)。
したがって `HasCoherentFunctional K` は**無条件に成り立つ**
(`ABC3.Found.PGC.hasCoherentFunctional`)。

★**有向集合上の全射逆系は極限が空になりうる**ので、可算共終列はどうしても要る
(本ファイルの §3 は ℕ 上の従属選択である)。
-/

namespace ABC3.Found.PGC

namespace OrbitModel

/-! ## 1. ★★★抽象核 —— 有限部分群上の軌道和(語彙ゼロ)

`avgSum N a = ∑_{n ∈ N} n • a`。分岐・付値・Galois どころか、体すら出てこない。 -/

section AvgSum

variable {G : Type*} [Group G] {A : Type*} [AddCommGroup A] [DistribMulAction G A]

/-- **有限部分群 `N` 上の軌道和** `∑_{n ∈ N} n • a`。 -/
def avgSum (N : Subgroup G) [Fintype N] (a : A) : A := ∑ n : N, (n : G) • a

theorem avgSum_add (N : Subgroup G) [Fintype N] (a b : A) :
    avgSum N (a + b) = avgSum N a + avgSum N b := by
  simp [avgSum, smul_add, Finset.sum_add_distrib]

@[simp] theorem avgSum_zero (N : Subgroup G) [Fintype N] : avgSum N (0 : A) = 0 := by
  simp [avgSum]

/-- 軌道和は `N` で不変。 -/
theorem avgSum_smul_mem (N : Subgroup G) [Fintype N] (a : A) (m : G) (hm : m ∈ N) :
    m • avgSum N a = avgSum N a := by
  show m • (∑ n : N, (n : G) • a) = ∑ n : N, (n : G) • a
  rw [Finset.smul_sum]
  refine Fintype.sum_bijective (fun n : N => (⟨m, hm⟩ : N) * n) (Group.mulLeft_bijective _) _ _ ?_
  intro n
  show m • ((n : G) • a) = ((⟨m, hm⟩ * n : N) : G) • a
  rw [smul_smul]
  rfl

/-- ★`N` が正規なら軌道和は `G`-同変。 -/
theorem avgSum_smul (N : Subgroup G) [hN : N.Normal] [Fintype N] (a : A) (g : G) :
    g • avgSum N a = avgSum N (g • a) := by
  show g • (∑ n : N, (n : G) • a) = ∑ n : N, (n : G) • (g • a)
  rw [Finset.smul_sum]
  refine Fintype.sum_bijective
    (fun n : N => (⟨g * (n : G) * g⁻¹, hN.conj_mem (n : G) n.2 g⟩ : N)) ?_ _ _ ?_
  · rw [Function.bijective_iff_has_inverse]
    exact ⟨fun n => ⟨g⁻¹ * (n : G) * g, by simpa using hN.conj_mem (n : G) n.2 g⁻¹⟩,
      fun n => Subtype.ext (by show g⁻¹ * (g * (n : G) * g⁻¹) * g = (n : G); group),
      fun n => Subtype.ext (by show g * (g⁻¹ * (n : G) * g) * g⁻¹ = (n : G); group)⟩
  · intro n
    show g • ((n : G) • a) = (g * (n : G) * g⁻¹) • (g • a)
    rw [smul_smul, smul_smul]
    congr 1
    group

/-- `N` が固定する元の上では軌道和は `|N|` 倍。 -/
theorem avgSum_of_fixed (N : Subgroup G) [Fintype N] (a : A) (h : ∀ n ∈ N, n • a = a) :
    avgSum N a = (Fintype.card N) • a := by
  have hn : ∀ n : N, (n : G) • a = a := fun n => h (n : G) n.2
  show (∑ n : N, (n : G) • a) = _
  rw [Finset.sum_congr rfl (fun n _ => hn n)]
  simp [Finset.card_univ]

end AvgSum

/-! ## 2. ★★★抽象核 —— Maschke 上げ(語彙ゼロ) -/

section Maschke

variable {G : Type*} [Group G] {A : Type*} [AddCommGroup A] [DistribMulAction G A]
variable {M : Type*} [AddCommGroup M]
variable {Q : Type*} [Group Q] {B : Type*} [AddCommGroup B] [DistribMulAction Q B]

/-- ★**良い汎関数どうしは同変加法同型で移り合う**(Frobenius 相互律の直接の帰結)。

`Φ_{μ₀}` と `Φ_μ` がともに全単射なら `φ := Φ_{μ₀}⁻¹ ∘ Φ_μ` が求めるもの。 -/
theorem exists_addEquiv_of_bijective (mu0 mu : B →+ M)
    (h0 : Function.Bijective (orbitHom (G := Q) mu0))
    (h : Function.Bijective (orbitHom (G := Q) mu)) :
    ∃ phi : B ≃+ B, (∀ (q : Q) (b : B), phi (q • b) = q • phi b) ∧ ∀ b, mu0 (phi b) = mu b := by
  classical
  set E0 := AddEquiv.ofBijective (orbitHom (G := Q) mu0) h0 with hE0
  set E := AddEquiv.ofBijective (orbitHom (G := Q) mu) h with hE
  have hE0a : ∀ b : B, (E0 b : Q → M) = orbitHom (G := Q) mu0 b := fun _ => rfl
  have hEa : ∀ b : B, (E b : Q → M) = orbitHom (G := Q) mu b := fun _ => rfl
  have key : ∀ b : B, orbitHom (G := Q) mu0 (E0.symm (E b)) = orbitHom (G := Q) mu b := by
    intro b
    rw [← hE0a, E0.apply_symm_apply, hEa]
  refine ⟨E.trans E0.symm, ?_, ?_⟩
  · intro q b
    refine E0.injective ?_
    show orbitHom (G := Q) mu0 (E0.symm (E (q • b))) = orbitHom (G := Q) mu0 (q • E0.symm (E b))
    rw [key]
    funext x
    show orbitFun mu (q • b) x = orbitFun mu0 (q • E0.symm (E b)) x
    rw [orbitFun_smul, orbitFun_smul]
    exact congrFun (key b).symm (q⁻¹ * x)
  · intro b
    have hb := congrFun (key b) 1
    simpa using hb

/-- ★★★**Maschke 上げ(抽象核)**。

`ι : B →+ A` が全射 `π : G ↠ Q` に沿って同変で、**同変な収縮** `ρ : A →+ B`
(`ρ ∘ ι = id`)を持つとする。このとき `A` 上の良い汎関数 `λ` と `B` 上の良い汎関数 `μ` が
あれば、`λ` を取り替えて **`μ` を延長する良い汎関数**が取れる。

★証明は初等的で、`Ψ := ι ∘ φ ∘ ρ + (1 − ι ∘ ρ)` が同変自己同型であること
(逆は `φ⁻¹` 版)を見るだけである。★**群環も中心冪等元も単元の持ち上げも使わない。** -/
theorem exists_bijective_extension (pi : G →* Q) (iota : B →+ A) (rho : A →+ B)
    (hiota : ∀ (g : G) (b : B), iota (pi g • b) = g • iota b)
    (hrho : ∀ (g : G) (a : A), rho (g • a) = pi g • rho a)
    (hri : ∀ b : B, rho (iota b) = b)
    (lam : A →+ M) (hlam : Function.Bijective (orbitHom (G := G) lam))
    (hlev : Function.Bijective (orbitHom (G := Q) (lam.comp iota)))
    (mu : B →+ M) (hmu : Function.Bijective (orbitHom (G := Q) mu)) :
    ∃ lam' : A →+ M, Function.Bijective (orbitHom (G := G) lam') ∧
      ∀ b : B, lam' (iota b) = mu b := by
  classical
  obtain ⟨phi, hphieq, hphival⟩ := exists_addEquiv_of_bijective (lam.comp iota) mu hlev hmu
  set corr : (B ≃+ B) → A → A := fun f a => iota (f (rho a)) + a - iota (rho a) with hcorr
  have hrho_corr : ∀ (f : B ≃+ B) (a : A), rho (corr f a) = f (rho a) := by
    intro f a
    simp [hcorr, map_sub, map_add, hri]
  have hcorr_corr : ∀ (f : B ≃+ B) (a : A), corr f (corr f.symm a) = a := by
    intro f a
    have h1 : rho (corr f.symm a) = f.symm (rho a) := hrho_corr f.symm a
    show iota (f (rho (corr f.symm a))) + corr f.symm a - iota (rho (corr f.symm a)) = a
    rw [h1]
    show iota (f (f.symm (rho a))) + (iota (f.symm (rho a)) + a - iota (rho a))
        - iota (f.symm (rho a)) = a
    rw [f.apply_symm_apply]
    abel
  have hcorr_add : ∀ (f : B ≃+ B) (x y : A), corr f (x + y) = corr f x + corr f y := by
    intro f x y
    show iota (f (rho (x + y))) + (x + y) - iota (rho (x + y))
        = (iota (f (rho x)) + x - iota (rho x)) + (iota (f (rho y)) + y - iota (rho y))
    rw [map_add, map_add, map_add, map_add]
    abel
  set Psi : A →+ A := AddMonoidHom.mk' (corr phi) (hcorr_add phi) with hPsi
  have hPsibij : Function.Bijective Psi := by
    constructor
    · intro x y hxy
      have hx := hcorr_corr phi.symm x
      have hy := hcorr_corr phi.symm y
      rw [phi.symm_symm] at hx hy
      rw [← hx, ← hy]
      show corr phi.symm (corr phi x) = corr phi.symm (corr phi y)
      exact congrArg _ hxy
    · intro a
      exact ⟨corr phi.symm a, hcorr_corr phi a⟩
  have hPsismul : ∀ (g : G) (a : A), Psi (g • a) = g • Psi a := by
    intro g a
    show iota (phi (rho (g • a))) + g • a - iota (rho (g • a))
        = g • (iota (phi (rho a)) + a - iota (rho a))
    rw [hrho g a, hphieq, hiota, hiota, smul_sub, smul_add]
  have hPsiiota : ∀ b : B, Psi (iota b) = iota (phi b) := by
    intro b
    show iota (phi (rho (iota b))) + iota b - iota (rho (iota b)) = iota (phi b)
    rw [hri]
    abel
  refine ⟨lam.comp Psi, ?_, ?_⟩
  · have hcompeq : (orbitHom (G := G) (lam.comp Psi) : A → (G → M))
        = (orbitHom (G := G) lam) ∘ Psi := by
      funext a
      funext g
      show lam (Psi (g⁻¹ • a)) = lam (g⁻¹ • Psi a)
      rw [hPsismul]
    rw [show (orbitHom (G := G) (lam.comp Psi) : A → (G → M)) = _ from hcompeq]
    exact hlam.comp hPsibij
  · intro b
    show lam (Psi (iota b)) = mu b
    rw [hPsiiota]
    exact hphival b

end Maschke

/-! ## 3. ★★★抽象核 —— ℕ 上の従属選択と直極限(語彙ゼロ) -/

section Tower

/-- ★**ℕ 上の従属選択**。各段で 1 歩進めるなら、全段の族が取れる。 -/
theorem exists_seq_of_step {X : ℕ → Type*} (R : ∀ n, X n → X (n + 1) → Prop)
    (x0 : X 0) (hstep : ∀ (n : ℕ) (x : X n), ∃ y, R n x y) :
    ∃ x : (∀ n, X n), ∀ n, R n (x n) (x (n + 1)) := by
  classical
  choose g hg using hstep
  refine ⟨fun n => Nat.rec (motive := X) x0 g n, ?_⟩
  intro n
  exact hg n _

/-- ★**直極限**。増大する部分加法群の列で覆われ、段ごとに整合する汎関数の族があれば、
全体の汎関数が 1 本取れる。 -/
theorem exists_addMonoidHom_of_tower {A M : Type*} [AddCommGroup A] [AddCommGroup M]
    (S : ℕ → AddSubgroup A) (hmono : Monotone S) (hcov : ∀ a : A, ∃ n, a ∈ S n)
    (f : ∀ n, S n →+ M)
    (hcompat : ∀ (n : ℕ) (a : A) (h : a ∈ S n) (h' : a ∈ S (n + 1)),
      f (n + 1) ⟨a, h'⟩ = f n ⟨a, h⟩) :
    ∃ g : A →+ M, ∀ (n : ℕ) (a : A) (h : a ∈ S n), g a = f n ⟨a, h⟩ := by
  classical
  set F : ℕ → A → M := fun n a => if h : a ∈ S n then f n ⟨a, h⟩ else 0 with hF
  have hFval : ∀ (n : ℕ) (a : A) (h : a ∈ S n), F n a = f n ⟨a, h⟩ := by
    intro n a h
    simp [hF, h]
  have hFadd : ∀ (n : ℕ) (x y : A), x ∈ S n → y ∈ S n →
      F n (x + y) = F n x + F n y := by
    intro n x y hx hy
    rw [hFval n x hx, hFval n y hy, hFval n (x + y) (add_mem hx hy)]
    exact map_add (f n) ⟨x, hx⟩ ⟨y, hy⟩
  have key : ∀ (m n : ℕ), n ≤ m → ∀ a : A, a ∈ S n → F m a = F n a := by
    intro m
    induction m with
    | zero =>
      intro n hn a _
      obtain rfl : n = 0 := Nat.le_zero.1 hn
      rfl
    | succ m ih =>
      intro n hn a ha
      rcases Nat.lt_or_ge m n with hlt | hge
      · obtain rfl : n = m + 1 := by omega
        rfl
      · have ham : a ∈ S m := hmono hge ha
        rw [hFval (m + 1) a (hmono (Nat.le_succ m) ham),
          hcompat m a ham (hmono (Nat.le_succ m) ham), ← hFval m a ham]
        exact ih n hge a ha
  have hcovN : ∀ a : A, a ∈ S (Nat.find (hcov a)) := fun a => Nat.find_spec (hcov a)
  refine ⟨AddMonoidHom.mk' (fun a => F (Nat.find (hcov a)) a) ?_, ?_⟩
  · intro x y
    have hx : F (max (max (Nat.find (hcov x)) (Nat.find (hcov y))) (Nat.find (hcov (x + y)))) x
        = F (Nat.find (hcov x)) x :=
      key _ _ (le_trans (le_max_left _ _) (le_max_left _ _)) x (hcovN x)
    have hy : F (max (max (Nat.find (hcov x)) (Nat.find (hcov y))) (Nat.find (hcov (x + y)))) y
        = F (Nat.find (hcov y)) y :=
      key _ _ (le_trans (le_max_right _ _) (le_max_left _ _)) y (hcovN y)
    have hxy : F (max (max (Nat.find (hcov x)) (Nat.find (hcov y))) (Nat.find (hcov (x + y))))
        (x + y) = F (Nat.find (hcov (x + y))) (x + y) :=
      key _ _ (le_max_right _ _) _ (hcovN _)
    show F (Nat.find (hcov (x + y))) (x + y) = F (Nat.find (hcov x)) x + F (Nat.find (hcov y)) y
    rw [← hxy, ← hx, ← hy]
    exact hFadd _ x y
      (hmono (le_trans (le_max_left _ _) (le_max_left _ _)) (hcovN x))
      (hmono (le_trans (le_max_right _ _) (le_max_left _ _)) (hcovN y))
  · intro n a h
    show F (Nat.find (hcov a)) a = f n ⟨a, h⟩
    rw [← hFval n a h]
    have h1 : F (max n (Nat.find (hcov a))) a = F n a := key _ _ (le_max_left _ _) a h
    have h2 : F (max n (Nat.find (hcov a))) a = F (Nat.find (hcov a)) a :=
      key _ _ (le_max_right _ _) a (hcovN a)
    rw [← h2, h1]

end Tower

end OrbitModel

/-! ## 4. 具体層 —— 2 つの有限水準の間の配線

★**設計**: 中間体を 2 層またぐ `Algebra ↥L ↥L'` を**作らない**
(`lean-idioms.md` #59:2 層をまたぐ `rfl` は kernel を止める)。
代わりに**すべての射を `Γ_K` を経由して定義する**——
`Gal(L'/K) ↠ Gal(L/K)` は `Γ_K ↠ Gal(L'/K)` の余核性(`MonoidHom.liftOfSurjective`)で作る。 -/

section TwoLevels

open ABC3.Skeleton.PGC OrbitModel

variable {p : ℕ} [Fact p.Prime]

attribute [local instance] carrier_charZero

variable (K : PAdicLocalField p) (L L' : IntermediateField K.carrier K.closure)

/-- 水準の包含 `↥L →+ ↥L'`。 -/
def levelIncl (h : L ≤ L') : L →+ L' where
  toFun b := ⟨(b : K.closure), h b.2⟩
  map_zero' := rfl
  map_add' _ _ := rfl

@[simp] theorem coe_levelIncl (h : L ≤ L') (b : L) :
    ((levelIncl K L L' h b : L') : K.closure) = (b : K.closure) := rfl

theorem levelIncl_injective (h : L ≤ L') : Function.Injective (levelIncl K L L' h) :=
  fun _ _ hxy => Subtype.ext (congrArg (fun z : L' => (z : K.closure)) hxy)

/-- 包含は `levelIota` と両立する(★定義から `rfl`)。 -/
theorem levelIota_comp_levelIncl (h : L ≤ L') :
    (levelIota K L').comp (levelIncl K L L' h) = levelIota K L := rfl

variable [Normal K.carrier L] [Normal K.carrier L']

/-- ★水準の間の制限射 `Gal(L'/K) ↠ Gal(L/K)`。★`Γ_K` を経由して作る。 -/
noncomputable def levelRes (h : L ≤ L') : (L' ≃ₐ[K.carrier] L') →* (L ≃ₐ[K.carrier] L) :=
  (levelPi K L').liftOfSurjective (levelPi_surjective K L') ⟨levelPi K L, by
    intro g hg
    rw [MonoidHom.mem_ker] at hg ⊢
    rw [levelPi_eq_one_iff] at hg ⊢
    rw [IntermediateField.mem_fixingSubgroup_iff] at hg ⊢
    exact fun x hx => hg x (h hx)⟩

@[simp] theorem levelRes_levelPi (h : L ≤ L') (g : K.absGal) :
    levelRes K L L' h (levelPi K L' g) = levelPi K L g :=
  MonoidHom.liftOfRightInverse_comp_apply _ _ _ _ g

theorem levelRes_surjective (h : L ≤ L') : Function.Surjective (levelRes K L L' h) := by
  intro σ
  obtain ⟨g, rfl⟩ := levelPi_surjective K L σ
  exact ⟨levelPi K L' g, levelRes_levelPi K L L' h g⟩

/-- ★遷移の同変性(2 水準版)。 -/
theorem levelIncl_smul (h : L ≤ L') (g : L' ≃ₐ[K.carrier] L') (b : L) :
    levelIncl K L L' h (levelRes K L L' h g • b) = g • levelIncl K L L' h b := by
  obtain ⟨γ, rfl⟩ := levelPi_surjective K L' g
  apply Subtype.ext
  rw [levelRes_levelPi]
  show ((levelPi K L γ b : L) : K.closure)
      = ((levelPi K L' γ (levelIncl K L L' h b) : L') : K.closure)
  rw [coe_levelPi_apply, coe_levelPi_apply, coe_levelIncl]

/-- 制限射の核は下の水準を各点で固定する。 -/
theorem levelIncl_fixed (h : L ≤ L') (n : L' ≃ₐ[K.carrier] L')
    (hn : levelRes K L L' h n = 1) (b : L) :
    n • levelIncl K L L' h b = levelIncl K L L' h b := by
  obtain ⟨γ, rfl⟩ := levelPi_surjective K L' n
  rw [levelRes_levelPi] at hn
  have hγ : γ ∈ L.fixingSubgroup := (levelPi_eq_one_iff K L γ).1 hn
  apply Subtype.ext
  show ((levelPi K L' γ (levelIncl K L L' h b) : L') : K.closure) = (b : K.closure)
  rw [coe_levelPi_apply, coe_levelIncl]
  exact (IntermediateField.mem_fixingSubgroup_iff L γ).1 hγ (b : K.closure) b.2

/-- 制限射の核で固定される元は下の水準に属する(Galois 降下、2 水準版)。 -/
theorem mem_range_levelIncl_of_fixed (h : L ≤ L') (a : L')
    (hfix : ∀ n : L' ≃ₐ[K.carrier] L', levelRes K L L' h n = 1 → n • a = a) :
    a ∈ Set.range (levelIncl K L L' h) := by
  have ha : (a : K.closure) ∈ IntermediateField.fixedField L.fixingSubgroup := by
    rw [IntermediateField.mem_fixedField_iff]
    intro γ hγ
    have h1 : levelRes K L L' h (levelPi K L' γ) = 1 := by
      rw [levelRes_levelPi, levelPi_eq_one_iff]
      exact hγ
    have h2 := congrArg (fun z : L' => (z : K.closure)) (hfix _ h1)
    show γ (a : K.closure) = (a : K.closure)
    rw [← coe_levelPi_apply K L' γ a]
    exact h2
  rw [InfiniteGalois.fixedField_fixingSubgroup] at ha
  exact ⟨⟨(a : K.closure), ha⟩, Subtype.ext rfl⟩

/-- スカラーと `Γ_K` の作用は交換する。 -/
theorem absGal_smul_comm (γ : K.absGal) (c : K.carrier) (y : K.closure) :
    γ (c • y) = c • γ y := by
  rw [Algebra.smul_def, Algebra.smul_def, map_mul, AlgEquiv.commutes]

/-- ★★★**具体層の Maschke 上げ**:下の水準の良い汎関数は、上の水準の良い汎関数に延長できる。 -/
theorem exists_good_extension (h : L ≤ L')
    [FiniteDimensional K.carrier L] [FiniteDimensional K.carrier L'] [IsGalois K.carrier L']
    (mu : L →+ K.carrier)
    (hmu : Function.Bijective (orbitHom (G := L ≃ₐ[K.carrier] L) mu)) :
    ∃ lam : L' →+ K.carrier,
      Function.Bijective (orbitHom (G := L' ≃ₐ[K.carrier] L') lam) ∧
      ∀ b : L, lam (levelIncl K L L' h b) = mu b := by
  classical
  haveI : Fintype (L' ≃ₐ[K.carrier] L') := Fintype.ofFinite _
  haveI : Fintype ((levelRes K L L' h).ker) := Fintype.ofFinite _
  obtain ⟨lam0, hlam0⟩ := exists_levelwise_functional K L'
  have hmemL : ∀ x : L', ((avgSum (levelRes K L L' h).ker x : L') : K.closure) ∈ L := by
    intro x
    obtain ⟨b, hb⟩ := mem_range_levelIncl_of_fixed K L L' h (avgSum (levelRes K L L' h).ker x)
      (fun n hn => avgSum_smul_mem _ x n (MonoidHom.mem_ker.2 hn))
    rw [← hb]
    exact b.2
  set c : K.carrier := ((Fintype.card ((levelRes K L L' h).ker) : ℕ) : K.carrier) with hc
  have hc0 : c ≠ 0 := by
    rw [hc]
    exact Nat.cast_ne_zero.2 Fintype.card_ne_zero
  obtain ⟨rho, coe_rho⟩ : ∃ rho : L' →+ L, ∀ x : L',
      ((rho x : L) : K.closure) = c⁻¹ • ((avgSum (levelRes K L L' h).ker x : L') : K.closure) := by
    refine ⟨AddMonoidHom.mk'
      (fun x => c⁻¹ • (⟨((avgSum (levelRes K L L' h).ker x : L') : K.closure), hmemL x⟩ : L))
      ?_, fun _ => rfl⟩
    intro x y
    have hxy : (⟨((avgSum (levelRes K L L' h).ker (x + y) : L') : K.closure), hmemL (x + y)⟩ : L)
        = ⟨((avgSum (levelRes K L L' h).ker x : L') : K.closure), hmemL x⟩
          + ⟨((avgSum (levelRes K L L' h).ker y : L') : K.closure), hmemL y⟩ := by
      apply Subtype.ext
      show ((avgSum (levelRes K L L' h).ker (x + y) : L') : K.closure)
          = ((avgSum (levelRes K L L' h).ker x : L') : K.closure)
            + ((avgSum (levelRes K L L' h).ker y : L') : K.closure)
      rw [avgSum_add]
      rfl
    rw [hxy, smul_add]
  have hri : ∀ b : L, rho (levelIncl K L L' h b) = b := by
    intro b
    apply Subtype.ext
    rw [coe_rho]
    have h1 : avgSum (levelRes K L L' h).ker (levelIncl K L L' h b)
        = (Fintype.card ((levelRes K L L' h).ker)) • levelIncl K L L' h b :=
      avgSum_of_fixed _ _ (fun n hn => levelIncl_fixed K L L' h n (MonoidHom.mem_ker.1 hn) b)
    rw [h1]
    have h2 : (((Fintype.card ((levelRes K L L' h).ker)) • levelIncl K L L' h b : L')
        : K.closure) = (Fintype.card ((levelRes K L L' h).ker)) • ((b : K.closure)) := by
      simp
    rw [h2, ← Nat.cast_smul_eq_nsmul K.carrier, ← hc, inv_smul_smul₀ hc0]
  have hrhosmul : ∀ (g : L' ≃ₐ[K.carrier] L') (x : L'),
      rho (g • x) = levelRes K L L' h g • rho x := by
    intro g x
    obtain ⟨γ, rfl⟩ := levelPi_surjective K L' g
    apply Subtype.ext
    rw [coe_rho, levelRes_levelPi]
    show c⁻¹ • ((avgSum (levelRes K L L' h).ker (levelPi K L' γ • x) : L') : K.closure)
        = ((levelPi K L γ (rho x) : L) : K.closure)
    rw [coe_levelPi_apply, coe_rho, absGal_smul_comm, ← avgSum_smul]
    congr 1
    exact coe_levelPi_apply K L' γ (avgSum (levelRes K L L' h).ker x)
  have hlev : Function.Bijective
      (orbitHom (G := L ≃ₐ[K.carrier] L) (lam0.comp (levelIncl K L L' h))) :=
    bijective_descent (levelRes K L L' h) (levelRes_surjective K L L' h) (levelIncl K L L' h)
      (levelIncl_injective K L L' h) (levelIncl_smul K L L' h)
      (mem_range_levelIncl_of_fixed K L L' h) lam0 hlam0
  exact exists_bijective_extension (levelRes K L L' h) (levelIncl K L L' h) rho
    (levelIncl_smul K L L' h) hrhosmul hri lam0 hlam0 hlev mu hmu

end TwoLevels

/-! ## 5. 具体層 —— 塔の組み立て -/

section Assemble

open ABC3.Skeleton.PGC OrbitModel

variable {p : ℕ} [Fact p.Prime]

attribute [local instance] carrier_charZero

variable (K : PAdicLocalField p)

/-- ★★**可算共終 Galois 塔**:`K̄/K` の有限次 Galois 部分拡大の増大列で `K̄` を覆うもの。

★本ファイルが `HasCoherentFunctional K` に残す**唯一**の穴。
★**証明は `Found/PGC/AlgebraicGenerators.lean::hasCofinalGaloisTower`**(Krasner の補題)。 -/
def HasCofinalGaloisTower (K : PAdicLocalField p) : Prop :=
  ∃ T : ℕ → IntermediateField K.carrier K.closure, Monotone T ∧
    (∀ n, FiniteDimensional K.carrier (T n)) ∧ (∀ n, IsGalois K.carrier (T n)) ∧
    (∀ x : K.closure, ∃ n, x ∈ T n)

/-- 有限次部分拡大は塔のどこかに収まる(有限生成性から)。 -/
theorem le_of_cofinal_tower (T : ℕ → IntermediateField K.carrier K.closure) (hmono : Monotone T)
    (hcov : ∀ x : K.closure, ∃ n, x ∈ T n) (L : IntermediateField K.carrier K.closure)
    [FiniteDimensional K.carrier L] : ∃ n, L ≤ T n := by
  classical
  have hfg : L.FG := IntermediateField.essFiniteType_iff.1 inferInstance
  obtain ⟨t, htfin, hadj⟩ := IntermediateField.fg_def.1 hfg
  haveI : Finite t := htfin.to_subtype
  choose ν hν using fun x : t => hcov (x : K.closure)
  obtain ⟨n, hn⟩ := Finite.exists_le ν
  refine ⟨n, ?_⟩
  rw [← hadj, IntermediateField.adjoin_le_iff]
  intro x hx
  exact hmono (hn ⟨x, hx⟩) (hν ⟨x, hx⟩)

/-- ★★★**塔から両立汎関数**。 -/
theorem hasCoherentFunctional_of_cofinalTower (h : HasCofinalGaloisTower K) :
    HasCoherentFunctional K := by
  classical
  obtain ⟨T, hmono, hfd, hgal, hcov⟩ := h
  haveI : ∀ n, FiniteDimensional K.carrier (T n) := hfd
  haveI : ∀ n, IsGalois K.carrier (T n) := hgal
  obtain ⟨mu0, hmu0⟩ := exists_levelwise_functional K (T 0)
  obtain ⟨mu, hmu⟩ := OrbitModel.exists_seq_of_step
    (X := fun n => {f : (T n) →+ K.carrier //
      Function.Bijective (orbitHom (G := T n ≃ₐ[K.carrier] T n) f)})
    (R := fun n x y => ∀ b : (T n),
      y.1 (levelIncl K (T n) (T (n + 1)) (hmono (Nat.le_succ n)) b) = x.1 b)
    ⟨mu0, hmu0⟩
    (by
      intro n x
      obtain ⟨lam, hlam, hval⟩ := exists_good_extension K (T n) (T (n + 1))
        (hmono (Nat.le_succ n)) x.1 x.2
      exact ⟨⟨lam, hlam⟩, hval⟩)
  obtain ⟨g, hg⟩ := OrbitModel.exists_addMonoidHom_of_tower
    (S := fun n => (T n).toSubfield.toAddSubgroup)
    (fun n m hnm => fun x hx => hmono hnm hx)
    (fun a => hcov a)
    (fun n => (mu n).1)
    (by
      intro n a ha ha'
      exact hmu n ⟨a, ha⟩)
  refine ⟨g, ?_⟩
  intro L hfdL hgalL
  haveI := hfdL
  haveI := hgalL
  obtain ⟨n, hn⟩ := le_of_cofinal_tower K T hmono hcov L
  have hTn : Function.Bijective
      (orbitHom (G := T n ≃ₐ[K.carrier] T n) (g.comp (levelIota K (T n)))) := by
    have heq : g.comp (levelIota K (T n)) = (mu n).1 := by
      ext b
      exact hg n (b : K.closure) b.2
    rw [heq]
    exact (mu n).2
  exact bijective_descent (levelRes K L (T n) hn) (levelRes_surjective K L (T n) hn)
    (levelIncl K L (T n) hn) (levelIncl_injective K L (T n) hn) (levelIncl_smul K L (T n) hn)
    (mem_range_levelIncl_of_fixed K L (T n) hn) (g.comp (levelIota K (T n))) hTn

/-! ### 5.1 ★★塔の十分条件 —— 「可算個の生成元」

`HasCofinalGaloisTower` は可算性の主張で、その**古典的な出所**は
「`K̄` の各元は、ある可算集合の元 1 つで `K` 上生成される」である
(Krasner の補題:`f` に十分近い `g` の根は同じ体を生成する)。
★実際にそれが成り立つことは `Found/PGC/AlgebraicGenerators.lean` で証明した
(`hasCountableGenerators`。★係数を `K` の**可算稠密部分集合**から取る)。
★歴史的な形として「ℚ 上代数的な元で生成される」版も残しておく。 -/

/-- ★★**代数的生成元**:`K̄` の各元は、**ℚ 上代数的**な元 1 つで `K` 上生成される部分体に入る。

★古典的には Krasner の補題の帰結(「`f` に十分近い `g` の根は同じ体を生成する」+
「ℚ は ℚ_p で稠密」)。★本木では、これより弱い
`HasCountableGenerators`(可算稠密集合から係数を取る)で十分だったので、
そちらを `AlgebraicGenerators.lean` で証明した。★この形は**証明していない**
(ℚ の稠密性を `K` まで持ち上げる手間が余分にかかるため)。 -/
def HasAlgebraicGenerators (K : PAdicLocalField p) : Prop :=
  ∀ x : K.closure, ∃ α : K.closure,
    IsAlgebraic ℚ α ∧ x ∈ IntermediateField.adjoin K.carrier {α}

/-- ★★**可算個の生成元**:`K̄` の各元が、ある**可算集合**の元 1 つで `K` 上生成される。 -/
def HasCountableGenerators (K : PAdicLocalField p) : Prop :=
  ∃ D : Set K.closure, D.Countable ∧
    ∀ x : K.closure, ∃ α ∈ D, x ∈ IntermediateField.adjoin K.carrier {α}

/-- ℚ 上代数的な元は可算個しかない(`Algebraic.countable`)。 -/
theorem hasCountableGenerators_of_algebraicGenerators (h : HasAlgebraicGenerators K) :
    HasCountableGenerators K :=
  ⟨{x : K.closure | IsAlgebraic ℚ x}, Algebraic.countable (R := ℚ) (A := K.closure),
    fun x => by
      obtain ⟨α, hα, hx⟩ := h x
      exact ⟨α, hα, hx⟩⟩

/-- ★★★**可算個の生成元から可算共終塔**。

生成元を `f : ℕ → K̄` で数え上げ、`U (n+1) := U n ⊔ (f (n+1) の Galois 閉包)` と積み上げる。 -/
theorem hasCofinalGaloisTower_of_countableGenerators (h : HasCountableGenerators K) :
    HasCofinalGaloisTower K := by
  classical
  obtain ⟨D, hcount, hgen⟩ := h
  have hne : D.Nonempty := by
    obtain ⟨α, hα, -⟩ := hgen 0
    exact ⟨α, hα⟩
  obtain ⟨f, hf⟩ := hcount.exists_eq_range hne
  obtain ⟨U, hU0, hUsucc⟩ : ∃ U : ℕ → FiniteGaloisIntermediateField K.carrier K.closure,
      U 0 = FiniteGaloisIntermediateField.adjoin K.carrier {f 0} ∧
      ∀ k, U (k + 1) = U k ⊔ FiniteGaloisIntermediateField.adjoin K.carrier {f (k + 1)} :=
    ⟨fun n => Nat.rec (FiniteGaloisIntermediateField.adjoin K.carrier {f 0})
      (fun k acc => acc ⊔ FiniteGaloisIntermediateField.adjoin K.carrier {f (k + 1)}) n,
     rfl, fun _ => rfl⟩
  have hmono : Monotone U := by
    refine monotone_nat_of_le_succ (fun n => ?_)
    rw [hUsucc n]
    exact le_sup_left
  have hmemU : ∀ n, f n ∈ (U n).toIntermediateField := by
    intro n
    cases n with
    | zero =>
      rw [hU0]
      exact FiniteGaloisIntermediateField.subset_adjoin K.carrier {f 0} rfl
    | succ k =>
      have h1 : f (k + 1)
          ∈ (FiniteGaloisIntermediateField.adjoin K.carrier {f (k + 1)}).toIntermediateField :=
        FiniteGaloisIntermediateField.subset_adjoin K.carrier {f (k + 1)} rfl
      have h2 : FiniteGaloisIntermediateField.adjoin K.carrier {f (k + 1)} ≤ U (k + 1) := by
        rw [hUsucc k]
        exact le_sup_right
      exact (FiniteGaloisIntermediateField.le_iff _ _).1 h2 h1
  refine ⟨fun n => (U n).toIntermediateField, ?_, fun n => inferInstance, fun n => inferInstance,
    ?_⟩
  · intro n m hnm
    exact (FiniteGaloisIntermediateField.le_iff _ _).1 (hmono hnm)
  · intro x
    obtain ⟨α, hα, hx⟩ := hgen x
    have hmem : α ∈ Set.range f := by rw [← hf]; exact hα
    obtain ⟨n, rfl⟩ := hmem
    exact ⟨n, IntermediateField.adjoin_simple_le_iff.2 (hmemU n) hx⟩

/-- ★代数的生成元から可算共終塔(上の 2 本の合成)。 -/
theorem hasCofinalGaloisTower_of_algebraicGenerators (h : HasAlgebraicGenerators K) :
    HasCofinalGaloisTower K :=
  hasCofinalGaloisTower_of_countableGenerators K (hasCountableGenerators_of_algebraicGenerators K h)

/-- ★★★**[pGC] Proposition 2.1 への配線**(塔の仮説つき)。 -/
theorem prop_2_1_of_cofinalTower
    (h : ∀ K : PAdicLocalField p, HasCofinalGaloisTower K) :
    RecoverableAsAddModule (p := p) (fun K => K.closure) :=
  prop_2_1_of_hasCoherentFunctional
    (fun K => hasCoherentFunctional_of_cofinalTower K (h K))

end Assemble

end ABC3.Found.PGC
