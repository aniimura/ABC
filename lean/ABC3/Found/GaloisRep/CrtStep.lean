import ABC3.Found.GaloisRep.PrimePowerEquiv

/-!
# Galois (G1) 第 70 ブロック —— **★★★★★中国剰余定理の段**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★一般の `n` を素数冪へ落とす

`n = m·k`(`gcd(m,k) = 1`)、`n•x = 0` のとき

    A[m] × A[k]  ≃+  A,   (x, y) ↦ x + y

★単射: `x + y = 0` なら `x = −y ∈ A[m] ∩ A[k]`、Bezout `um + vk = 1` から `x = 0`。
★★全射: `z = (vk)•z + (um)•z` で、`(vk)•z ∈ A[m]`、`(um)•z ∈ A[k]`。

## ★★部分群の中の捩れ

`d ∣ N` なら `A[N][d] = A[d]` である(`d•x = 0 ⟹ N•x = 0`)。
★これで素数冪の段(第 69)の仮定が `A` の仮定から引き継げる。

## ★★★`ZMod` 側

    (ℤ/m)² × (ℤ/k)²  ≃+  (ℤ/mk)²

★`ZMod.chineseRemainder` と `AddEquiv.prodProdProdComm` を繋ぐだけ。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `addEquiv_crt` | ★★★★★**`A[m] × A[k] ≃+ A`** |
| `kerSubEquiv` | ★★`A[N][d] ≃+ A[d]`(`d ∣ N`) |
| `zmodProdCrt` | ★★`(ℤ/m)² × (ℤ/k)² ≃+ (ℤ/mk)²` |
-/

namespace ABC3.Found.GaloisRep

universe u

/-- ★★★★★**互いに素な分解** `A[m] × A[k] ≃+ A`。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to -/
theorem addEquiv_crt {A : Type u} [AddCommGroup A] (m k : ℕ) (hco : Nat.Coprime m k)
    (hexp : ∀ x : A, (m * k) • x = 0) :
    Nonempty (((nsmulHom A m).ker × (nsmulHom A k).ker) ≃+ A) := by
  obtain ⟨u, v, huv⟩ : ∃ u v : ℤ, u * (m : ℤ) + v * (k : ℤ) = 1 := by
    have h1 : Int.gcd (m : ℤ) (k : ℤ) = 1 := by simpa [Int.gcd_natCast_natCast] using hco
    exact ⟨Int.gcdA (m : ℤ) (k : ℤ), Int.gcdB (m : ℤ) (k : ℤ), by
      have hg := Int.gcd_eq_gcd_ab (m : ℤ) (k : ℤ)
      rw [h1] at hg
      push_cast at hg
      linarith [hg]⟩
  set f : ((nsmulHom A m).ker × (nsmulHom A k).ker) →+ A :=
    ((nsmulHom A m).ker.subtype.comp (AddMonoidHom.fst _ _))
      + ((nsmulHom A k).ker.subtype.comp (AddMonoidHom.snd _ _)) with hf
  have hfval : ∀ z : ((nsmulHom A m).ker × (nsmulHom A k).ker), f z = (z.1 : A) + (z.2 : A) := by
    intro z; rw [hf]; simp
  have hinj : Function.Injective f := by
    rw [injective_iff_map_eq_zero]
    rintro ⟨⟨x, hx⟩, ⟨y, hy⟩⟩ hz
    rw [hfval] at hz
    simp only at hz
    have hxm : (m : ℤ) • x = 0 := by rw [natCast_zsmul]; exact hx
    have hyk : (k : ℤ) • y = 0 := by rw [natCast_zsmul]; exact hy
    have hxy : x = -y := by linear_combination (norm := abel) hz
    have hxk : (k : ℤ) • x = 0 := by rw [hxy, zsmul_neg, hyk, neg_zero]
    have hx0 : x = 0 := by
      have hh : ((u * (m : ℤ) + v * (k : ℤ)) : ℤ) • x = x := by rw [huv, one_zsmul]
      rw [add_zsmul, mul_zsmul, mul_zsmul, hxm, hxk, zsmul_zero, zsmul_zero, add_zero] at hh
      exact hh.symm
    have hy0 : y = 0 := by rw [hxy] at hx0; simpa using hx0
    ext <;> simp [hx0, hy0]
  have hsurj : Function.Surjective f := by
    intro z
    have hx : (m : ℕ) • ((v * (k : ℤ)) • z) = 0 := by
      rw [← natCast_zsmul ((v * (k : ℤ)) • z) m, smul_smul,
        show (m : ℤ) * (v * (k : ℤ)) = v * (((m * k : ℕ) : ℤ)) by push_cast; ring,
        mul_zsmul, natCast_zsmul, hexp z, zsmul_zero]
    have hy : (k : ℕ) • ((u * (m : ℤ)) • z) = 0 := by
      rw [← natCast_zsmul ((u * (m : ℤ)) • z) k, smul_smul,
        show (k : ℤ) * (u * (m : ℤ)) = u * (((m * k : ℕ) : ℤ)) by push_cast; ring,
        mul_zsmul, natCast_zsmul, hexp z, zsmul_zero]
    refine ⟨⟨⟨(v * (k : ℤ)) • z, hx⟩, ⟨(u * (m : ℤ)) • z, hy⟩⟩, ?_⟩
    rw [hfval]
    simp only
    rw [← add_zsmul, show v * (k : ℤ) + u * (m : ℤ) = 1 by linarith [huv], one_zsmul]
  exact ⟨AddEquiv.ofBijective f ⟨hinj, hsurj⟩⟩

/-- ★★**部分群の中の `d` 捩れは元の `d` 捩れと同じ**(`d ∣ N`)。 -/
def kerSubEquiv {A : Type u} [AddCommGroup A] (N d : ℕ) (hd : d ∣ N) :
    (nsmulHom ((nsmulHom A N).ker) d).ker ≃+ (nsmulHom A d).ker where
  toFun := fun x => ⟨(x.1 : A), by
    have hx : d • (x.1 : (nsmulHom A N).ker) = 0 := x.2
    exact congrArg Subtype.val hx⟩
  invFun := fun y => ⟨⟨(y.1 : A), by
      obtain ⟨t, ht⟩ := hd
      have hy : d • (y.1 : A) = 0 := y.2
      show N • (y.1 : A) = 0
      rw [ht, mul_comm, mul_smul, hy, smul_zero]⟩, by
    show d • (⟨(y.1 : A), _⟩ : (nsmulHom A N).ker) = 0
    ext
    exact y.2⟩
  left_inv := by intro x; ext; rfl
  right_inv := by intro y; ext; rfl
  map_add' := by intro x y; rfl

/-- ★★**`(ℤ/m)² × (ℤ/k)² ≃+ (ℤ/mk)²`**(互いに素)。 -/
def zmodProdCrt (m k : ℕ) (hco : Nat.Coprime m k) :
    ((ZMod m × ZMod m) × (ZMod k × ZMod k)) ≃+ (ZMod (m * k) × ZMod (m * k)) :=
  (AddEquiv.prodProdProdComm (ZMod m) (ZMod m) (ZMod k) (ZMod k)).trans
    (((ZMod.chineseRemainder hco).toAddEquiv.symm).prodCongr
      ((ZMod.chineseRemainder hco).toAddEquiv.symm))

/-! ## ★出典の紐付け(`.src`) -/

def addEquiv_crt.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(捩れ部分群の中国剰余分解)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
