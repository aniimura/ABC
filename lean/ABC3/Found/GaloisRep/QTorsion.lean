import ABC3.Found.GaloisRep.TateLimit
import Mathlib.RingTheory.RootsOfUnity.AlgebraicallyClosed
import Mathlib.FieldTheory.IsAlgClosed.Basic
import Mathlib.GroupTheory.QuotientGroup.Basic

/-!
# Galois (G6) 第 103 ブロック —— **★★★★★★`Kˣ/q^ℤ` の捩れと Tate 加群**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★Tate 一意化の**右辺**を先に固める

`Skeleton/GaloisRep/TateUniformization.lean` が待っているのは

    E_q(K) ≅ Kˣ / q^ℤ

である。★その**左辺**(級数 `X(u,q)`, `Y(u,q)` と加法公式)は葉 (a)–(e) として残っているが、
★★**右辺は今すぐ全部やれる**——`Kˣ/q^ℤ` は純粋に群論の対象だからである。

## ★★★★本ブロックが出す等式

`ζ` を原始 `N` 乗根、`π^N = q`、`q` が 1 の冪根でないとすると

    (Kˣ/q^ℤ)[N]  ≃+  (ℤ/N)²        ——  生成元は `ζ` と `π`

★単射性は「`ζ^m π^k ∈ q^ℤ` ⟹ `N ∣ m`, `N ∣ k`」であり、
`q` が 1 の冪根でないことだけを使う。
★★全射性は「`x^N ∈ q^ℤ` ⟹ `x·π^{-j}` が 1 の `N` 乗根」であり、
mathlib の `IsPrimitiveRoot.zpowers_eq`(体の `N` 乗根は `ζ` で生成)を使う。

## ★★★★★★そして Tate 加群が出る

★第 73–77 ブロックの `addEquiv_limTors` は
「すべての `m ≥ 1` で `A[m] ≃ (ℤ/m)²`」から `T_l A ≃+ ℤ_l²` を出す装置であった。
★★代数閉・標数 0 なら `π` も `ζ` も**常に存在する**ので、その仮定が満たされる:

    T_l (Kˣ/q^ℤ)  ≃+  ℤ_l × ℤ_l

★★★★★★**これは Tate 曲線の Tate 加群そのものである**
——一意化 `E_q(K) ≅ Kˣ/q^ℤ` が繋がった瞬間に (G5) の像計算へ渡る。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `qTorsHom` | ★`(a,b) ↦ ζ^a·π^b` |
| `qTorsHom_injective` | ★★★単射(`q` が 1 の冪根でないこと**だけ**を使う) |
| `qTorsHom_ker_le_range` | ★★★全射(`N` 乗根が `ζ` で生成されること) |
| `qQuot_torsion_addEquiv` | ★★★★★**`(Kˣ/q^ℤ)[N] ≃+ (ℤ/N)²`** |
| `zpow_eq_zero_of_val` | ★★付値が自明でなければ 1 の冪根でない |
| `tateModule_qQuot` | ★★★★★★**`T_l(Kˣ/q^ℤ) ≃+ ℤ_l²`** |
-/

namespace ABC3.Found.GaloisRep

open QuotientGroup

/-! ## ★`ZMod N` からの持ち上げ -/

/-- ★`g^N ∈ q^ℤ` なら `a ↦ g^a` は `ZMod N → Kˣ/q^ℤ` を定める。 -/
noncomputable def qLift {K : Type} [Field K] (N : ℕ) (q g : Kˣ)
    (hg : g ^ N ∈ Subgroup.zpowers q) :
    ZMod N →+ Additive (Kˣ ⧸ Subgroup.zpowers q) :=
  ZMod.lift N ⟨zmultiplesHom _ (Additive.ofMul (QuotientGroup.mk' (Subgroup.zpowers q) g)), by
    simp only [zmultiplesHom_apply, ← ofMul_zpow, zpow_natCast, ← map_pow]
    rw [ofMul_eq_zero]
    exact (QuotientGroup.eq_one_iff _).2 hg⟩

theorem qLift_intCast {K : Type} [Field K] (N : ℕ) (q g : Kˣ)
    (hg : g ^ N ∈ Subgroup.zpowers q) (m : ℤ) :
    qLift N q g hg (m : ZMod N)
      = Additive.ofMul (QuotientGroup.mk (s := Subgroup.zpowers q) (g ^ m)) := by
  rw [qLift, ZMod.lift_coe]
  simp only [zmultiplesHom_apply, ← ofMul_zpow]
  rfl

/-! ## ★★`(a, b) ↦ ζ^a · π^b` -/

/-- ★★`(a, b) ↦ ζ^a · π^b` —— `(ℤ/N)² → Kˣ/q^ℤ`。 -/
noncomputable def qTorsHom {K : Type} [Field K] (N : ℕ) (q π ζ : Kˣ)
    (hπ : π ^ N ∈ Subgroup.zpowers q) (hζ : ζ ^ N ∈ Subgroup.zpowers q) :
    (ZMod N × ZMod N) →+ Additive (Kˣ ⧸ Subgroup.zpowers q) :=
  AddMonoidHom.mk' (fun x => qLift N q ζ hζ x.1 + qLift N q π hπ x.2) (by
    intro a b
    show qLift N q ζ hζ (a.1 + b.1) + qLift N q π hπ (a.2 + b.2) = _
    rw [map_add, map_add]
    abel)

theorem qTorsHom_intCast {K : Type} [Field K] (N : ℕ) (q π ζ : Kˣ)
    (hπ : π ^ N ∈ Subgroup.zpowers q) (hζ : ζ ^ N ∈ Subgroup.zpowers q) (m k : ℤ) :
    qTorsHom N q π ζ hπ hζ ((m : ZMod N), (k : ZMod N))
      = Additive.ofMul (QuotientGroup.mk (s := Subgroup.zpowers q) (ζ ^ m * π ^ k)) := by
  show qLift N q ζ hζ (m : ZMod N) + qLift N q π hπ (k : ZMod N) = _
  rw [qLift_intCast, qLift_intCast]
  rfl

/-! ## ★★★単射性 —— `q` が 1 の冪根でないことだけを使う -/

/-- ★★★**単射性**。`ζ^m·π^k ∈ q^ℤ` なら `N ∣ m` かつ `N ∣ k`。

★証明: `q = π^N` を使うと `π^(N j − k) = ζ^m`。★★両辺を `N` 乗すると
`ζ^N = 1` から `q^(N j − k) = 1`、★★★`q` が 1 の冪根でないので `k = N j`。
★★★★すると `ζ^m = π^0 = 1` となり `ζ` の位数から `N ∣ m`。 -/
theorem qTorsHom_injective {K : Type} [Field K] {N : ℕ} {q π ζ : Kˣ}
    {hπ : π ^ N ∈ Subgroup.zpowers q} {hζ : ζ ^ N ∈ Subgroup.zpowers q}
    (hπq : π ^ N = q) (hζ1 : ζ ^ N = 1)
    (hord : ∀ m : ℤ, ζ ^ m = 1 → (N : ℤ) ∣ m)
    (hqinf : ∀ j : ℤ, q ^ j = 1 → j = 0) :
    Function.Injective (qTorsHom N q π ζ hπ hζ) := by
  have hqp : ∀ s : ℤ, q ^ s = π ^ ((N : ℤ) * s) := by
    intro s
    rw [zpow_mul, zpow_natCast, hπq]
  rw [injective_iff_map_eq_zero]
  rintro ⟨a, b⟩ h
  obtain ⟨m, rfl⟩ := ZMod.intCast_surjective a
  obtain ⟨k, rfl⟩ := ZMod.intCast_surjective b
  rw [show ((m : ZMod N), (k : ZMod N)) = ((m : ZMod N), (k : ZMod N)) from rfl,
    qTorsHom_intCast, ofMul_eq_zero, QuotientGroup.eq_one_iff] at h
  obtain ⟨j, hj⟩ := h
  have hj' : q ^ j = ζ ^ m * π ^ k := hj
  have h1 : π ^ ((N : ℤ) * j - k) = ζ ^ m := by
    rw [zpow_sub, ← hqp, hj', mul_inv_cancel_right]
  have h2 : q ^ ((N : ℤ) * j - k) = 1 := by
    have hr := congrArg (fun x : Kˣ => x ^ (N : ℤ)) h1
    simp only [← zpow_mul] at hr
    rw [hqp, mul_comm ((N : ℤ)) (((N : ℤ) * j - k)), hr, mul_comm, zpow_mul, zpow_natCast,
      hζ1, one_zpow]
  have h3 : (N : ℤ) * j - k = 0 := hqinf _ h2
  have hk : k = (N : ℤ) * j := by omega
  have hm : (N : ℤ) ∣ m := by
    apply hord
    rw [← h1, h3, zpow_zero]
  refine Prod.ext ?_ ?_
  · simpa using (ZMod.intCast_zmod_eq_zero_iff_dvd m N).2 hm
  · have hz : ((k : ℤ) : ZMod N) = 0 := by
      rw [(ZMod.intCast_zmod_eq_zero_iff_dvd k N)]
      exact ⟨j, hk⟩
    simpa using hz

/-! ## ★★★像は `N` 捩れ全体 -/

/-- ★像は `N` 捩れに入る(`(ℤ/N)²` では `N • x = 0` だから)。 -/
theorem qTorsHom_mem_ker {K : Type} [Field K] {N : ℕ} {q π ζ : Kˣ}
    {hπ : π ^ N ∈ Subgroup.zpowers q} {hζ : ζ ^ N ∈ Subgroup.zpowers q}
    (x : ZMod N × ZMod N) :
    qTorsHom N q π ζ hπ hζ x
      ∈ (nsmulHom (Additive (Kˣ ⧸ Subgroup.zpowers q)) N).ker := by
  rw [mem_ker_nsmulHom, ← map_nsmul]
  have hx : N • x = 0 := by
    refine Prod.ext ?_ ?_ <;> simp
  rw [hx, map_zero]

/-- ★★★**全射性**。`x^N ∈ q^ℤ` なら `x·π^{-j}` は 1 の `N` 乗根であり、
`ζ` の冪である。 -/
theorem qTorsHom_ker_le_range {K : Type} [Field K] {N : ℕ} {q π ζ : Kˣ}
    {hπ : π ^ N ∈ Subgroup.zpowers q} {hζ : ζ ^ N ∈ Subgroup.zpowers q}
    (hπq : π ^ N = q)
    (hgen : ∀ u : Kˣ, u ^ N = 1 → u ∈ Subgroup.zpowers ζ)
    (y : Additive (Kˣ ⧸ Subgroup.zpowers q))
    (hy : y ∈ (nsmulHom (Additive (Kˣ ⧸ Subgroup.zpowers q)) N).ker) :
    ∃ x, qTorsHom N q π ζ hπ hζ x = y := by
  have hqp : ∀ s : ℤ, q ^ s = π ^ ((N : ℤ) * s) := by
    intro s
    rw [zpow_mul, zpow_natCast, hπq]
  rw [mem_ker_nsmulHom] at hy
  obtain ⟨x, hx⟩ := QuotientGroup.mk_surjective (s := Subgroup.zpowers q) (Additive.toMul y)
  have hy1 : ((Additive.toMul y) ^ N) = 1 := by
    rw [← ofMul_eq_zero, ofMul_pow]
    simpa using hy
  have hx1 : ((x ^ N : Kˣ) : Kˣ ⧸ Subgroup.zpowers q) = 1 := by
    rw [← QuotientGroup.mk'_apply, map_pow, QuotientGroup.mk'_apply, hx]
    exact hy1
  obtain ⟨j, hj⟩ := (QuotientGroup.eq_one_iff _).1 hx1
  have hj' : q ^ j = x ^ N := hj
  have hxN : x ^ N = π ^ ((N : ℤ) * j) := by rw [← hqp, hj']
  have hexp : (x * (π ^ j)⁻¹) ^ N = x ^ N * ((π ^ j) ^ N)⁻¹ := by
    rw [mul_pow, inv_pow]
  have huN : (x * (π ^ j)⁻¹) ^ N = 1 := by
    rw [hexp, hxN, ← zpow_natCast (π ^ j) N, ← zpow_mul, mul_comm j (N : ℤ), mul_inv_cancel]
  obtain ⟨i, hi⟩ := hgen _ huN
  have hi' : ζ ^ i = x * (π ^ j)⁻¹ := hi
  refine ⟨((i : ZMod N), (j : ZMod N)), ?_⟩
  rw [qTorsHom_intCast]
  have hxi : ζ ^ i * π ^ j = x := by
    rw [hi', inv_mul_cancel_right]
  rw [hxi]
  show Additive.ofMul (QuotientGroup.mk x) = y
  rw [hx]
  rfl

/-! ## ★★★★★`(Kˣ/q^ℤ)[N] ≃+ (ℤ/N)²` -/

/-- ★★★★★**`(Kˣ/q^ℤ)[N] ≃+ (ℤ/N)²`**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem qQuot_torsion_addEquiv {K : Type} [Field K] {N : ℕ} {q π ζ : Kˣ}
    (hπq : π ^ N = q) (hζ1 : ζ ^ N = 1)
    (hord : ∀ m : ℤ, ζ ^ m = 1 → (N : ℤ) ∣ m)
    (hgen : ∀ u : Kˣ, u ^ N = 1 → u ∈ Subgroup.zpowers ζ)
    (hqinf : ∀ j : ℤ, q ^ j = 1 → j = 0) :
    Nonempty ((ZMod N × ZMod N)
      ≃+ (nsmulHom (Additive (Kˣ ⧸ Subgroup.zpowers q)) N).ker) := by
  have hπ : π ^ N ∈ Subgroup.zpowers q := by rw [hπq]; exact Subgroup.mem_zpowers q
  have hζ : ζ ^ N ∈ Subgroup.zpowers q := by rw [hζ1]; exact one_mem _
  refine ⟨AddEquiv.ofBijective
    ((qTorsHom N q π ζ hπ hζ).codRestrict _ (fun x => qTorsHom_mem_ker x)) ⟨?_, ?_⟩⟩
  · intro a b hab
    exact qTorsHom_injective (hπ := hπ) (hζ := hζ) hπq hζ1 hord hqinf (congrArg Subtype.val hab)
  · rintro ⟨y, hy⟩
    obtain ⟨x, hxy⟩ := qTorsHom_ker_le_range (hπ := hπ) (hζ := hζ) hπq hgen y hy
    exact ⟨x, Subtype.ext hxy⟩

/-- ★★★★**原始 `N` 乗根版**——mathlib の `IsPrimitiveRoot` に接続する。 -/
theorem qQuot_torsion_addEquiv_of_primitiveRoot {K : Type} [Field K] {N : ℕ} [NeZero N]
    {q π ζ : Kˣ} (hπq : π ^ N = q) (hζ : IsPrimitiveRoot ζ N)
    (hqinf : ∀ j : ℤ, q ^ j = 1 → j = 0) :
    Nonempty ((ZMod N × ZMod N)
      ≃+ (nsmulHom (Additive (Kˣ ⧸ Subgroup.zpowers q)) N).ker) := by
  refine qQuot_torsion_addEquiv hπq hζ.pow_eq_one
    (fun m hm => (hζ.zpow_eq_one_iff_dvd m).1 hm) (fun u hu => ?_) hqinf
  rw [hζ.zpowers_eq]
  exact (mem_rootsOfUnity N u).2 hu

/-! ## ★★「`q` は 1 の冪根でない」を付値から出す -/

/-- ★★**付値が自明でないなら `q` は 1 の冪根ではない**。

★Tate 曲線の設定(`q` が離散付値環の極大イデアルに入り `q ≠ 0`)では
`v(q) > 0` なので、この補題が仮定を供給する。 -/
theorem zpow_eq_zero_of_val {K : Type} [Field K] (v : Kˣ →* Multiplicative ℤ) (q : Kˣ)
    (hv : v q ≠ 1) (j : ℤ) (hj : q ^ j = 1) : j = 0 := by
  have h := congrArg v hj
  rw [map_zpow, map_one] at h
  have h2 : j * Multiplicative.toAdd (v q) = 0 := by
    have hc := congrArg Multiplicative.toAdd h
    rwa [toAdd_zpow, toAdd_one, smul_eq_mul] at hc
  rcases mul_eq_zero.1 h2 with h4 | h4
  · exact h4
  · exact absurd (toAdd_eq_zero.mp h4) hv

/-! ## ★★★★★★Tate 加群 -/

/-- ★★★★代数閉・標数 0 なら、`q` が 1 の冪根でない限り
`(Kˣ/q^ℤ)[m] ≃+ (ℤ/m)²` が**すべての `m ≥ 1`** で成り立つ。

★`π` は `IsAlgClosed.exists_pow_nat_eq` で、
★★`ζ` は `HasEnoughRootsOfUnity`(分離閉体の instance)で得る。 -/
theorem qQuot_torsion_card {K : Type} [Field K] [IsAlgClosed K] [CharZero K] {q : Kˣ}
    (hqinf : ∀ j : ℤ, q ^ j = 1 → j = 0) (m : ℕ) (hm : 1 ≤ m) :
    Nonempty ((ZMod m × ZMod m)
      ≃+ (nsmulHom (Additive (Kˣ ⧸ Subgroup.zpowers q)) m).ker) := by
  haveI : NeZero m := ⟨by omega⟩
  haveI : NeZero ((m : ℕ) : K) := ⟨Nat.cast_ne_zero.2 (by omega)⟩
  obtain ⟨z, hz⟩ := IsAlgClosed.exists_pow_nat_eq (q : K) (n := m) (by omega)
  have hz0 : z ≠ 0 := by
    intro h
    rw [h, zero_pow (by omega)] at hz
    exact q.ne_zero hz.symm
  have hπq : (Units.mk0 z hz0) ^ m = q := by
    refine Units.ext ?_
    push_cast
    exact hz
  obtain ⟨ζ0, hζ0⟩ := HasEnoughRootsOfUnity.exists_primitiveRoot K m
  exact qQuot_torsion_addEquiv_of_primitiveRoot hπq (hζ0.isUnit_unit (by omega)) hqinf

/-- ★★★★★★**`T_l(Kˣ/q^ℤ) ≃+ ℤ_l²`** —— Tate 曲線の Tate 加群。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★第 73–77 ブロックの `addEquiv_limTors` に `qQuot_torsion_card` を流し込むだけである。
★★★一意化 `E_q(K) ≅ Kˣ/q^ℤ`(`Skeleton/GaloisRep/TateUniformization.lean`)が
繋がった瞬間、これが Tate 曲線の Tate 加群になる。 -/
theorem tateModule_qQuot {K : Type} [Field K] [IsAlgClosed K] [CharZero K] {q : Kˣ}
    (hqinf : ∀ j : ℤ, q ^ j = 1 → j = 0) (l : ℕ) [Fact l.Prime] :
    Nonempty (limTors (Additive (Kˣ ⧸ Subgroup.zpowers q)) l ≃+ (ℤ_[l] × ℤ_[l])) := by
  refine addEquiv_limTors ?_ ?_ l
  · intro m hm
    haveI : NeZero m := ⟨by omega⟩
    obtain ⟨e⟩ := qQuot_torsion_card hqinf m hm
    exact Finite.of_equiv _ e.toEquiv
  · intro m hm
    obtain ⟨e⟩ := qQuot_torsion_card hqinf m hm
    rw [← Nat.card_congr e.toEquiv, Nat.card_prod, Nat.card_zmod, sq]

/-! ## ★出典の紐付け(`.src`) -/

def qQuot_torsion_addEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Kˣ/q^ℤ の N 捩れが (ℤ/N)² であること)",
    sectionId := "genell-def-3-3" }

def tateModule_qQuot.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Kˣ/q^ℤ の Tate 加群が ℤ_l² であること)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
