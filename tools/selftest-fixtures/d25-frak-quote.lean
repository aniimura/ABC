-- D25 fixture: フラクトゥール記法つきの正しい引用(通るべき)
--   IUTchIII は書体だけで別対象を区別する(MOD 対 𝔪𝔬𝔡、LGP 対 𝔩𝔤𝔭、log 対 𝔩𝔬𝔤)。
--   平文表現 X[frak] が照合射影で落ちることを確かめる。
namespace Fixture
/-- 原文 (IUTchIII p.154):
> (Ind1) the indeterminacies induced by the automorphisms of the procession
> of D[scr]^-prime-strips Prc(n,◦D[frak]^_T);
-/
theorem frakQuotedOk : True := trivial
def frakQuotedOk.src : ABC3.Meta.Source :=
  { paper := "IUTchIII", pdfPage := 154, item := "Theorem 3.11, (i), (Ind1)",
    sectionId := "thm-3-11-ind1" }
def frakQuotedOk.needs : List ABC3.Meta.ProofObligation := []
end Fixture
