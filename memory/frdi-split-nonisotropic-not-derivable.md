---
name: frdi-split-nonisotropic-not-derivable
description: [FrdI] Prop 2.5 (iii) の「分裂は非 isotropic でも使える」は Def 2.3 (a)(b) からは導けなかった。単元部分の全射性が欠けている。
metadata:
  type: project
---

[FrdI] `Proposition 2.5, (iii)` の証明(物理 p.49)で原文は、非 isotropic な `A` でも
特性分裂 `𝒪^×(A) × τ(A) ≅ 𝒪^▷(A)` が使えると述べ、根拠に
**`Definition 2.3, (a), (b)` だけ**を挙げる。★**2026-08-17 に機械で追い、導けなかった。**

**詰まる場所**: `x ∈ 𝒪^▷(A)` を単射 `ι : 𝒪^▷(A) ↪ 𝒪^▷(A^istr)`(`Prop 2.2, (iv)`)で
送ると (a) で `ι x = u' · t'` と分裂し、(b) で `t' = ι t` なる `t` が一意に取れる
(`hullTau_existsUnique`、`lean/ABC3/Found/FrdI/Prop25.lean`)。
★**しかし `u' ∈ 𝒪^×(A^istr)` が `ι` の像に入ることが出ない。**
`Prop 2.2, (iv)` は「モノイドの自然な単射」しか主張せず、単元部分の全射性は言わない。
`ι t` は一般に可逆でないので `u' = ι x · (ι t)⁻¹` とも書けない。

**Why:** これは `Proposition 2.5` を数えられるかどうかの**律速**である。この穴があるため
`Prop25iii.lean` の `psiMap` は `IsOfIsotropicType P` を仮定したままで、
得られているのは実質「`𝒞^istr` 上の (iii)」である(圏同値そのものは通っている)。
★原文の段 8(isotropic hull が単射 ⟹ 一般の `𝒞` へ)も、**Ψ を非 isotropic な対象で
定義する**ために同じ分裂を要求するので、迂回にならない。

**How to apply:**
- ★**まだ ③ `sourceGap` ではない**。③ を主張するには「(a)(b) を満たすが非 isotropic な
  対象で分裂が破れる Frobenioid」という **falsifier** が要る。現状は「導けなかった」測定。
- 埋める候補: (α) `𝒪^×(A) → 𝒪^×(A^istr)` の全射性を別経路で出す
  (isotropic hull が mono であることと `Definition 1.3, (v), (a)` あたり)、
  (β) `Definition 2.3` に条件を足す —— ただし (β) は既に数えた項目の忠実度を下げるので、
  やるなら「原文が (a)(b) から従うとする内容を、我々は導けないので条件として明示する」と
  docstring に書くこと。
- 関連: [[frdi-prop21-quantifier-false]](原文の literal な主張が偽だった別の例)。
